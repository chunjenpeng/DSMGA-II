/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <map>
#include <unordered_set>

#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"

#include <iomanip>
using namespace std;

DSMGA2::DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff) {


    previousFitnessMean = -INF;
    ell = n_ell;
    nCurrent = (n_nInitial/2)*2;  // has to be even

    Chromosome::function = (Chromosome::Function)fffff;
    Chromosome::nfe = 0;
    Chromosome::lsnfe = 0;
    Chromosome::hitnfe = 0;
    Chromosome::hit = false;

    selectionPressure = 2;
    maxGen = n_maxGen;
    maxFe = n_maxFe;

    graph.init(ell);

    bestIndex = -1;
    masks = new list<int>[ell];
    selectionIndex = new int[nCurrent];
    orderN = new int[nCurrent];
    orderELL = new int[ell];
    population = new Chromosome[nCurrent];
    fastCounting = new FastCounting[ell];

    for (int i = 0; i < ell; i++)
        fastCounting[i].init(nCurrent);


    pHash.clear();
    for (int i=0; i<nCurrent; ++i) {
        population[i].initR(ell);
        double f = population[i].getFitness();
        pHash[population[i].getKey()] = f;
    }

    if (GHC) {
        for (int i=0; i < nCurrent; i++)
            population[i].GHC();
    }
}


DSMGA2::~DSMGA2 () {
    delete []masks;
    delete []orderN;
    delete []orderELL;
    delete []selectionIndex;
    delete []population;
    delete []fastCounting;
}



bool DSMGA2::isSteadyState () {

    if (stFitness.getNumber () <= 0)
        return false;

    if (previousFitnessMean < stFitness.getMean ()) {
        previousFitnessMean = stFitness.getMean () + EPSILON;
        return false;
    }

    return true;
}



int DSMGA2::doIt (bool output) {
    generation = 0;
    lastMax = lastMean = lastMin = -INF;
    convergeCount = 0;
    RM_failed = RM_succeed = BM_failed = BM_succeed = 0;

    while (!shouldTerminate ()) {
        oneRun (output);
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {

    if (CACHE)
        Chromosome::cache.clear();
#ifdef DEBUG
    printPopulation();
#endif
    mixing();


    double max = -INF;
    stFitness.reset ();

    for (int i = 0; i < nCurrent; ++i) {
        double fitness = population[i].getFitness();
        if (fitness > max) {
            max = fitness;
            bestIndex = i;
        }
        stFitness.record (fitness);

    }

    if (output)
        showStatistics ();

    ++generation;
}

bool DSMGA2::converged() {
    if (stFitness.getMax() == lastMax &&
        stFitness.getMean() == lastMean &&
        stFitness.getMin() == lastMin)
        convergeCount++;
    else
        convergeCount = 0;

    lastMax = stFitness.getMax();
    lastMean = stFitness.getMean();
    lastMin = stFitness.getMin();
    
    return (convergeCount > 100) ? true : false;
}

bool DSMGA2::shouldTerminate () {
    bool termination = false;

    if (maxFe != -1) {
        if (Chromosome::nfe > maxFe)
            termination = true;
    }

    if (maxGen != -1) {
        if (generation > maxGen)
            termination = true;
    }


    if (population[0].getMaxFitness() <= stFitness.getMax() )
        termination = true;


    if (stFitness.getMax() - EPSILON <= stFitness.getMean() )
        termination = true;

    if (converged() )
        termination = true;

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {
    /*
#ifndef DEUBG 
    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
#endif    
#ifdef DEBUG
    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f nfe:%d\n ",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin (), Chromosome::nfe);
#endif
    */
    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f nfe:%d\n",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin (), Chromosome::nfe);
    printf("lsnfe:%d, nfe:%d, RM failed:%d, RM success:%d, BM failed:%d, BM success:%d\n"
        ,Chromosome::lsnfe, Chromosome::nfe, RM_failed, RM_succeed, BM_failed, BM_succeed );

    fflush(NULL);
}



void DSMGA2::buildFastCounting() {

    if (SELECTION) {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[selectionIndex[i]].getVal(j));
            }

    } else {
        for (int i = 0; i < nCurrent; i++)
            for (int j = 0; j < ell; j++) {
                fastCounting[j].setVal(i, population[i].getVal(j));
            }
    }

}

int DSMGA2::countOne(int x) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;

        val = fastCounting[x].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


int DSMGA2::countXOR(int x, int y) const {

    int n = 0;

    for (int i=0; i<fastCounting[0].lengthLong; ++i) {
        unsigned long val = 0;


        val = fastCounting[x].gene[i];

        val ^= fastCounting[y].gene[i];

        n += myBD.countOne(val);
    }

    return n;
}


void DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    bool evaluated = trial.isEvaluated();

    if (trial.getFitness() > des.getFitness()) {

        countSucceed(mask, des, evaluated);

        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        return;
    }

    countFailed(mask, des, evaluated);
}

void DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    bool evaluated = trial.isEvaluated();

    if (trial.getFitness() > des.getFitness()) {

        countSucceed(mask, des, evaluated);

        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        des = trial;
        return;
    }

    //2016-10-21
    if (trial.getFitness() > des.getFitness() - EPSILON) {
    //if (trial.getFitness() >= des.getFitness()) {

        countSucceed(mask, des, evaluated);

        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        des = trial;
        return;
    }

    countFailed(mask, des, evaluated);
}

bool DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    bool taken = false;
    list<int> sMask;
    vector< pair< list<int>, double > > sortedMasks; 

    //2016-11-11 
    sortMasks( ch, mask, sortedMasks );
    for( pair< list<int>, double > p : sortedMasks ){
        sMask = p.first;
        
        Chromosome trial(ell);
        trial = ch;

        for (list<int>::iterator it = sMask.begin(); it != sMask.end(); ++it) {
            trial.flip(*it);
        }
    
        #ifdef DEBUG
        cout << " Try Mask:";
        printMaskScore(p);
        #endif
        if (isInP(trial)){
            #ifdef DEBUG
            printMask(sMask);
            cout << "isInP";
            #endif
            break;
        }
    
        //2016-10-21
        //if (trial.getFitness() > ch.getFitness()) {
        //if (trial.getFitness() >= ch.getFitness()) {
        ++RM_failed;
        if (trial.getFitness() > ch.getFitness() - EPSILON) {
            #ifdef DEBUG 
            for (auto it = sortedMasks.begin(); it != sortedMasks.end(); ++it )
                printMaskScore( *it );
                cout << "\nTaken Mask:";
                printMask(sMask);
                cout << endl;

                printf("%11.6f before : ", ch.getFitness());
                for(int i = 0; i < ch.getLength(); i++)
                    cout << ch.getVal(i);
                cout << endl;
    
                printf("%11.6f after  : ", trial.getFitness());
                for(int i = 0; i < trial.getLength(); i++)
                    cout << trial.getVal(i);
                cout << endl;
            #endif

                --RM_failed;
                ++RM_succeed;

            pHash.erase(ch.getKey());
            pHash[trial.getKey()] = trial.getFitness();

            taken = true;
            ch = trial;
        }

        if (taken) break;
    }

    if (sMask.size() != 0) {
        mask = sMask;
    }
    
    return taken;

}

double DSMGA2::BMestimation( const Chromosome& ch, const list<int>& mask ){
    //Count patterns in Population
    map<string, int> counter;
    for (int n = 0; n < nCurrent; ++n) {
        string pattern;
        for (const int& i : mask){
            pattern += to_string(population[n].getVal(i));
        }
        ++counter[pattern];
    }

    vector< pair<string,int> > mapcopy(counter.begin(), counter.end());
    sort( mapcopy.begin(), mapcopy.end(),
        [](const pair< string, int > &left, const pair< string, int > &right){
                return left.second < right.second;
        });

    //Count ratio of succeedPattern and originalPattern in population
    double occur = 0;
    string succeedPattern, originalPattern;
    for (const int& i : mask) { 
        int allel = ch.getVal(i);
        originalPattern += to_string(allel);
        succeedPattern += (allel == 0) ? "1" : "0";
    }

    bool passOriginalPattern = false, passSucceedPattern = false;
    for (const auto& it : mapcopy) {
        if (it.first == succeedPattern) {
            occur += it.second;
            passSucceedPattern = true;
        }
        else if (it.first == originalPattern) {
            occur += it.second;
            passOriginalPattern = true;
        }
        if ( passOriginalPattern && passSucceedPattern ) break;
    }
    return occur/nCurrent;
}

void DSMGA2::populationMaskStatus( const Chromosome& ch, const list<int>& mask ){
    
    for (const int& i : mask){
        int allel = (ch.getVal(i) == 1) ? 0 : 1;
        cout << allel;
    }
    cout << " -> ";
    for (const int& i : mask){
        cout << ch.getVal(i);
    }
    cout << endl;

    map<string, int> counter;
    for (int n = 0; n < nCurrent; ++n) {
        string pattern;
        for (const int& i : mask){
            pattern += to_string(population[n].getVal(i));
        }
        ++counter[pattern];
    }
    printMapOrder(counter);
}

bool DSMGA2::matchPattern(Chromosome& source, list<int>& mask, Chromosome& des) {
    for (auto const& i : mask) {
        if (source.getVal(i) != des.getVal(i))
            return false;
    }
    return true;
}

void DSMGA2::countSucceed(list<int>& mask, Chromosome& des, bool evaluated) {
    string pattern;
    for (const int& i : mask){
        pattern += to_string(des.getVal(i));
    }
    ++succeedPattern[pattern];
    if( !evaluated )
        ++BM_succeed;
}

void DSMGA2::countFailed(list<int>& mask, Chromosome& des, bool evaluated) {
    string pattern;
    for (const int& i : mask){
        pattern += to_string(des.getVal(i));
    }
    ++failedPattern[pattern];
    if( !evaluated )
        ++BM_failed;
}

void DSMGA2::printMapOrder(map<string, int>& m){
    vector< pair<string,int> > mapcopy(m.begin(), m.end());
    sort( mapcopy.begin(), mapcopy.end(),
        [](const pair< string, int > &left, const pair< string, int > &right){
                return left.second > right.second;
        });
    for (const auto& it : mapcopy)
        cout << it.first << ":" << it.second << endl;
}

void DSMGA2::restrictedMixing(Chromosome& ch) {
    //BM_failed = BM_succeed = 0;
    succeedPattern.clear();
    failedPattern.clear();

    int r = myRand.uniformInt(0, ell-1);
    list<int> mask = masks[r];
    
    size_t size = findSize(ch, mask);
    if (size > (size_t)ell/2)
        size = ell/2;
        
    // prune mask to exactly size
    while (mask.size() > size)
        mask.pop_back();

    if(mask.size() == 0) return;
    bool taken = restrictedMixing(ch, mask);
        
    EQ = true;
    if (taken) {
#ifdef DEBUG
        cout << "\nBefore BM:" << BMestimation(ch, mask) << endl;
        populationMaskStatus(ch, mask);
#endif
        genOrderN();

        for (int i=0; i<nCurrent; ++i) {

            if (EQ)
                backMixingE(ch, mask, population[orderN[i]]);
            else
                backMixing(ch, mask, population[orderN[i]]);
        }
#ifdef DEBUG
        cout << "\nsucceedPattern:" << endl;
        printMapOrder(succeedPattern);
        //for (auto const& it : succeedPattern) 
        //    cout << it.first << ":" << it.second << endl;

        cout << "\nfailedPattern:" << endl;
        printMapOrder(failedPattern);
        //for (auto const& it : failedPattern) 
        //    cout << it.first << ":" << it.second << endl;
        cout << "\nAfter BM:" << BMestimation(ch, mask) << endl;
        populationMaskStatus(ch, mask);

        printf("\nlsnfe:%d, nfe:%d, RM failed:%d, RM success:%d, BM failed:%d, BM success:%d\n"
            ,Chromosome::lsnfe, Chromosome::nfe, RM_failed, RM_succeed, BM_failed, BM_succeed );
        cin.sync();
        cin.get();
#endif
    }
}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask) const {

    DLLA candidate(nCurrent);
    for (int i=0; i<nCurrent; ++i)
        candidate.insert(i);

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

        int allele = ch.getVal(*it);

        for (DLLA::iterator it2 = candidate.begin(); it2 != candidate.end(); ++it2) {
            if (population[*it2].getVal(*it) == allele)
                candidate.erase(*it2);

            if (candidate.isEmpty())
                break;
        }

        if (candidate.isEmpty())
            break;

        ++size;
    }

    return size;


}

size_t DSMGA2::findSize(Chromosome& ch, list<int>& mask, Chromosome& ch2) const {

    size_t size = 0;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
        if (ch.getVal(*it) == ch2.getVal(*it)) break;
        ++size;
    }
    return size;
}

void DSMGA2::printPopulation(void) const {
    cout << "population:" << endl;
    for (int i = 0; i < nCurrent; ++i){
        cout << setw(16) << " ";
        for (int j = 0; j < ell; ++j){
            cout << population[i].getVal(j);
        }
        cout << endl;
    }
}

void DSMGA2::printMask( const list<int> &mask ){
    cout << "[";
    for( const auto &m : mask ){
        cout << m << "-";
    }
    cout << "\b]";
}

void DSMGA2::printMaskScore( const pair< list<int>, double > &p ){
    list<int> mask = p.first;
    double score = p.second;
    printf("%.6f : [", score);
    for( const auto &m : mask ){
        printf("%d-", m);
    }
    printf("\b]\n");
}

void DSMGA2::generateRestMask( const list<int> &mask, vector<int> &rest ){
    bool *inMask = new bool[ell]();
    for (auto it = mask.begin(); it != mask.end(); ++it) {
        inMask[*it] = true;
    }
    for (int i = 0; i < ell; ++i) {
        if (!inMask[i])
            rest.push_back(i);
    }
}

double DSMGA2::averageEdge( const list<int> &mask ) {
    vector<int> rest;
    generateRestMask( mask, rest );
    
    double intraCluster = 0.0;
    if (mask.size() > 1){
        for (auto it = mask.begin(); next(it,1) != mask.end(); ++it) {
            for (auto it1 = next(it,1); it1 != mask.end(); ++it1){
                intraCluster += graph( *it, *it1 );
                //cout << "(" << *it << "," <<  *it1 << ")";
            }
        }
        intraCluster /= (mask.size() * (mask.size()-1) / 2.0);
    }

    double restIntraCluster = 0.0;
    if (rest.size() > 1) {
        for (auto it = rest.begin(); next(it,1) != rest.end(); ++it) {
            for (auto it1 = next(it,1); it1 != rest.end(); ++it1){
                restIntraCluster += graph( *it, *it1 );
            }
        }
        restIntraCluster /= (rest.size() * (rest.size()-1) / 2.0);
    }
    
    double interCluster = 0.0;
    for (auto it = mask.begin(); it != mask.end(); ++it) {
        for (auto it1 = rest.begin(); it1 != rest.end(); ++it1){
            interCluster += graph( *it, *it1 );
        }
    }
    interCluster /= (mask.size() * rest.size());

    //double score = intraCluster - interCluster;
    double score = restIntraCluster + intraCluster - interCluster;
    return score;
}

double DSMGA2::DaviesBouldin_index( const list<int> &mask ) {
    //if (generation%2 == 0 ) return 1.0/mask.size();
    if (mask.size() == 1) return 1.0;
    vector<int> rest;
    generateRestMask( mask, rest );
    
    double *S = new double [ell](); // initialized all to 0
    for (auto it = mask.begin(); it != mask.end(); ++it) {
        for ( auto it1 = mask.begin(); it1 != mask.end(); ++it1 ) {
            if ( *it != *it1 ) 
                S[*it] += 1-graph( *it, *it1 );
                //cout << *it << "-" << *it1 << ":" << graph( *it, *it1 ) << endl;
            /*if ( generation == 2 ) {
                cout << *it << "-" << *it1 << ":" << graph( *it, *it1 ) << endl;
                cin.sync();
                cin.get();
            }*/
        }
        //cout << mask.size() << "(" << *it << "):" << S[*it] << "/";
        S[*it] = (fabs(S[*it]) < EPSILON) ? 0 : S[*it] /= (mask.size()-1);
        //cout << S[*it] << endl;
    }

    /*
    for (auto it = rest.begin(); it != rest.end(); ++it) {
        for ( auto it1 = rest.begin(); it1 != rest.end(); ++it1 ) {
            if ( *it != *it1 ) 
                S[*it] += graph( *it, *it1 );
        }
        S[*it] = ( fasb(S[*it]) < EPSILON) ? 0 : S[*it] /= rest.size();
    }
    */

    double DB = -INF;
    for (auto i = mask.begin(); i != mask.end(); ++i) {
        for (auto j = rest.begin(); j != rest.end(); ++j) {
            double Mij = (1-graph( *i, *j ));
            //double Rij = (S[*i] + S[*j])/Mij;
            if(Mij < 0) printf("M%d-%d:%f", *i, *j, Mij); 
            //double Rij = ( fabs(Mij - 1) < EPSILON ) ? log(S[*i]) - log(Mij) : log(S[*i]) + log(Mij);
            double Rij = ( fabs(Mij) < EPSILON ) ? S[*i] : S[*i]/Mij;
            DB = ( Rij > DB ) ? Rij : DB;
        }
    }
    delete []S;
    return DB;
}

double DSMGA2::Dunn_index( const list<int> &mask ) {
    double score = 0.0;
    return score;
}

double DSMGA2::silhouette_coefficient( const list<int> &mask ) {
    double score = 0.0;
    return score;
}

double DSMGA2::clusterScore( const list<int> &mask ) {
    //return averageEdge( mask );
    return 1.0/DaviesBouldin_index( mask );
    //return Dunn_index( mask );
    //return silhouette_coefficient( mask );
}

void DSMGA2::sortMasks( const Chromosome& ch, list<int>& mask, 
                        vector< pair< list<int>, double > >& sortedMasks) {

    sortedMasks.clear();

    while (mask.size() > 0) {
        double RMscore = clusterScore(mask);
        double BMscore = BMestimation( ch, mask );
        double score = pow(RMscore, BMscore);
#ifdef DEBUG
        printf("RM: %.6f, BM: %.6f, score = %8.6f", RMscore, BMscore, score);
        printMask(mask);
        cout << endl;
#endif

        sortedMasks.push_back( make_pair(mask, score) );
        mask.pop_back();
    }

    sort( sortedMasks.begin(), sortedMasks.end(),
        [](const pair< list<int>,double > &left, const pair< list<int>,double > &right){
            if (fabs(left.second-right.second)<EPSILON) // float equal comparison
                return left.first.size() < right.first.size(); // favor small mask when tie
            else
                return left.second > right.second; // favor small DB index
        });
    
}

void DSMGA2::mixing() {

    if (SELECTION)
        selection();

    //* really learn model
    buildFastCounting();
    buildGraph();

    for (int i=0; i<ell; ++i)
        findClique(i, masks[i]);

    int repeat = (ell>50)? ell/50: 1;

    for (int k=0; k<repeat; ++k) {

        genOrderN();
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
            if (Chromosome::hit) break;
        }
        if (Chromosome::hit) break;
    }

}

inline bool DSMGA2::isInP(const Chromosome& ch) const {

    unordered_map<unsigned long, double>::const_iterator it = pHash.find(ch.getKey());
    return (it != pHash.end());
}

inline void DSMGA2::genOrderN() {
    myRand.uniformArray(orderN, nCurrent, 0, nCurrent-1);
}

inline void DSMGA2::genOrderELL() {
    myRand.uniformArray(orderELL, ell, 0, ell-1);
}

void DSMGA2::buildGraph() {

    int *one = new int [ell];
    for (int i=0; i<ell; ++i) {
        one[i] = countOne(i);
    }

    for (int i=0; i<ell; ++i) {

        for (int j=i+1; j<ell; ++j) {

            int n00, n01, n10, n11;
            int nX =  countXOR(i, j);

            n11 = (one[i]+one[j]-nX)/2;
            n10 = one[i] - n11;
            n01 = one[j] - n11;
            n00 = nCurrent - n01 - n10 - n11;

            double p00 = (double)n00/(double)nCurrent;
            double p01 = (double)n01/(double)nCurrent;
            double p10 = (double)n10/(double)nCurrent;
            double p11 = (double)n11/(double)nCurrent;

            double linkage;
            linkage = computeMI(p00,p01,p10,p11);
            //printf(" MI: (%2d-%2d) : %f\n", i, j, linkage);
            //DMC
            //linkage = DMC(p00, p01, p10, p11);
            //printf("DMC: (%2d-%2d) : %f\n", i, j, linkage);
            graph.write(i,j,linkage);
        }
    }


    delete []one;

}

double DSMGA2::DMC (double p00, double p01, double p10, double p11) {
    double a00 = p00*p11;
    double a01 = p01*p10;
    return (a00 > a01) ? (a00 - a01) : (a01 - a00);
}

// from 1 to ell, pick by max edge
void DSMGA2::findClique(int startNode, list<int>& result) {

    result.clear();

    DLLA rest(ell);
    genOrderELL();
    for (int i=0; i<ell; ++i) {
        if (orderELL[i]==startNode)
            result.push_back(orderELL[i]);
        else
            rest.insert(orderELL[i]);
    }

    double *connection = new double[ell];

    for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
        connection[*iter] = graph(startNode, *iter);

    while (!rest.isEmpty()) {

        double max = -INF;
        int index = -1;
        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter) {
            if (max < connection[*iter]) {
                max = connection[*iter];
                index = *iter;
            }
        }

        rest.erase(index);
        result.push_back(index);

        for (DLLA::iterator iter = rest.begin(); iter != rest.end(); ++iter)
            connection[*iter] += graph(index, *iter);
    }


    delete []connection;

}

double DSMGA2::computeMI(double a00, double a01, double a10, double a11) const {

    double p0 = a00 + a01;
    double q0 = a00 + a10;
    double p1 = 1-p0;
    double q1 = 1-q0;

    double join = 0.0;
    if (a00 > EPSILON)
        join += a00*log(a00);
    if (a01 > EPSILON)
        join += a01*log(a01);
    if (a10 > EPSILON)
        join += a10*log(a10);
    if (a11 > EPSILON)
        join += a11*log(a11);

    double p = 0.0;
    if (p0 > EPSILON)
        p += p0*log(p0);
    if (p1 > EPSILON)
        p += p1*log(p1);


    double q = 0.0;
    if (q0 > EPSILON)
        q += q0*log(q0);
    if (q1 > EPSILON)
        q += q1*log(q1);

    return -p-q+join;

}


void DSMGA2::selection () {
    tournamentSelection ();
}


// tournamentSelection without replacement
void DSMGA2::tournamentSelection () {
    int i, j;

    int randArray[selectionPressure * nCurrent];

    for (i = 0; i < selectionPressure; i++)
        myRand.uniformArray (randArray + (i * nCurrent), nCurrent, 0, nCurrent - 1);

    for (i = 0; i < nCurrent; i++) {

        int winner = 0;
        double winnerFitness = -INF;

        for (j = 0; j < selectionPressure; j++) {
            int challenger = randArray[selectionPressure * i + j];
            double challengerFitness = population[challenger].getFitness();

            if (challengerFitness > winnerFitness) {
                winner = challenger;
                winnerFitness = challengerFitness;
            }

        }
        selectionIndex[i] = winner;
    }
}
