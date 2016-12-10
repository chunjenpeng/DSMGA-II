/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <map>

#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"

#include <iomanip>
using namespace std;
//#define DEBUG

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
    RM_succeed = RM_failed = BM_succeed = BM_failed = 0;

    while (!shouldTerminate ()) {
        oneRun (output);
    }
    return generation;
}


void DSMGA2::oneRun (bool output) {

    if (CACHE)
        Chromosome::cache.clear();

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


bool DSMGA2::shouldTerminate () {
    bool
    termination = false;

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

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {
    #ifdef DEBUG 
    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n ",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
    #endif
    #ifndef DEBUG 
    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f nfe:%d\n",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin (), Chromosome::nfe);
    printf("lsnfe:%d, nfe:%d, RM failed:%d, RM success:%d, BM failed:%d, BM success:%d\n"
        ,Chromosome::lsnfe, Chromosome::nfe, RM_failed, RM_succeed, BM_failed, BM_succeed );
    #endif
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

void DSMGA2::printPopulation() const {
    cout << "population:" << endl;
    for (int i = 0; i < nCurrent; ++i){
        //cout << setw(20) << " ";
        for (int j = 0; j < ell; ++j){
            cout << population[i].getVal(j);
        }
        cout << endl;
    }
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

void DSMGA2::countSucceed(list<int>& mask, Chromosome& des, bool evaluated) {
    string pattern;
    for (const int& i : mask){
        pattern += to_string(des.getVal(i));
    }
    ++succeedPattern[pattern];
    if (!evaluated)
        ++BM_succeed;
}

void DSMGA2::countFailed(list<int>& mask, Chromosome& des, bool evaluated) {
    string pattern;
    for (const int& i : mask){
        pattern += to_string(des.getVal(i));
    }
    ++failedPattern[pattern];
    if (!evaluated)
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


void DSMGA2::savePattern( const Chromosome& ch, list<int>& mask){
    map<int, int> pattern;
    for (const int& i : mask)
        pattern[i] = ch.getVal(i);

    //if ( !patternExists(pattern) )
        patternList.push_back( make_pair(0.0, pattern) );
}

void DSMGA2::printPattern(const map<int, int>& pattern) {
    for( int i = 0; i < ell; ++i ){
        if( pattern.find(i) != pattern.end() )
            cout << pattern.at(i);
        else
            cout << "."; 
    }
    cout << endl;
}
    
bool DSMGA2::contradictPattern( const map<int, int>& pattern1, 
                                const map<int, int>& pattern2 ) { 
    for ( const auto& pair : pattern1 ) {
        int pos = pair.first;
        int allel = pair.second;
        if( pattern2.find(pos) != pattern2.end() )
            if (pattern2.at(pos) != allel)
                return true;
    }
    return false;
}

bool DSMGA2::patternExists( const map<int, int>& pattern ) {
    for (const auto& score_pattern : patternList )
        if( pattern == score_pattern.second )
            return true;
    return false;
}

void DSMGA2::merge( map<int, int>& mergedPattern,
                    const map<int, int>& pattern ) { 
    for ( const auto& pair : pattern ) 
        mergedPattern[pair.first] = pair.second;
}

double DSMGA2::BMestimation( map<int, int>& pattern ) {
    map<string, int> counter;
    for (int n = 0; n < nCurrent; ++n) {
        string sPattern;
        for (const auto& p : pattern) {
            sPattern += to_string(population[n].getVal(p.first));
        }
        ++counter[sPattern];
    }

    vector< pair<string,int> > mapcopy(counter.begin(), counter.end());
    sort( mapcopy.begin(), mapcopy.end(),
        [](const pair< string, int > &left, const pair< string, int > &right){
                return left.second < right.second;
        });

    double occur = 0;
    string succeedPattern, originalPattern;
    for (const auto& p : pattern) {
        succeedPattern += to_string(p.second);
        originalPattern += (p.second == 0) ? "1" : "0";
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

    //for (const auto& it : mapcopy)
    //    cout << it.first << ":" << it.second << endl;
    //cout << succeedPattern << ":" << occur/nCurrent << endl;
    //cin.sync();
    //cin.get();

    return occur/nCurrent;
}

void DSMGA2::mergeMasks() {
    
    for (auto& score_pattern : patternList ) { 
       score_pattern.first = BMestimation(score_pattern.second); 
    }

    patternList.sort(
        [] (const pair<double, map<int, int> >& lhs, const pair<double, map<int, int> >& rhs){
            if (fabs(lhs.first-rhs.first) < EPSILON)
                return lhs.second.size() > rhs.second.size();
            else
                return lhs.first > rhs.first;
    });

    printf("Original Patterns:%lu\n", patternList.size());
    for (const auto& score_pattern : patternList) {
        printf("%.6f        : ", score_pattern.first);
        printPattern( score_pattern.second );
    }
    cin.sync();
    cin.get();

    for (auto it1 = patternList.begin(); next(it1,1) != patternList.end(); ++it1) {
        map<int, int>& pattern1 = it1->second;

        for (auto it2 = next(it1, 1); it2 != patternList.end(); ++it2) {
            map<int, int>& pattern2 = it2->second;
        
            if ( !contradictPattern( pattern1, pattern2 ) ) {
                map<int, int> newPattern( pattern1 );
                merge( newPattern, pattern2 );
                if ( !patternExists(newPattern) ) { 
                    double score = it1->first + it2->first;
                    //double score = BMestimation(newPattern);
                    patternList.push_back( make_pair(score, newPattern) );
                }
            }
        }
    }

    patternList.sort(
        [] (const pair<double, map<int, int> >& lhs, const pair<double, map<int, int> >& rhs){
            if (fabs(lhs.first-rhs.first) < EPSILON)
                return lhs.second.size() > rhs.second.size();
            else
                return lhs.first > rhs.first;
    });

    /*
    printf("New Patterns:%lu\n", patternList.size());
    for (const auto& score_pattern : patternList) {
        printf("%.6f : ", score_pattern.first);
        printPattern( score_pattern.second );
    }
    cin.sync();
    cin.get();
    */
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


    bool taken = restrictedMixing(ch, mask);
    if(taken) 
        savePattern( ch, mask );
/*
    EQ = true;
    if (taken) {
#ifdef DEBUG
        cout << "\nBefore BM:" << endl;
        populationMaskStatus(ch, mask);
#endif
        genOrderN();

        for (int i=0; i<nCurrent; ++i) {
            //if (!matchPattern(ch, mask, population[orderN[i]])) continue;

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
        cout << "\nAfter BM:" << endl;
        populationMaskStatus(ch, mask);

        //Chromosome::nfe = RM_failed + RM_succeed + BM_failed + BM_succeed + nCurrent
        printf("\nRM_succeed:%d, RM failed:%d, BM succeed:%d, BM failed:%d, nfe:%d, lsnfe:%d\n"
            ,RM_succeed, RM_failed, BM_succeed, BM_failed, Chromosome::nfe, Chromosome::lsnfe);
        cin.sync();
        cin.get();
#endif
    }
*/
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

void DSMGA2::backMixingO(Chromosome& des) {
    //int n = mergedPatterns.size();
    //int *orderO = new int[n];  
    //myRand.uniformArray(orderO, n, 0, n-1);

    for (auto it = patternList.begin(); it != patternList.end(); ++it) {
        double& score = it-> first;
        map<int, int>& pattern = it->second;
        //printf("%.6f    Try : ", score);
        //printPattern( pattern );

        Chromosome trial(ell);
        trial = des;

        for( const auto& pos_allel : pattern )
            trial.setVal( pos_allel.first, pos_allel.second );
         
        bool evaluated = trial.isEvaluated();

        //if (trial.getFitness() > des.getFitness()) {
        if (trial.getFitness() > des.getFitness() - EPSILON) {
            if (!evaluated) ++BM_succeed;
            //#ifdef DEBUG
            printf("%.6f   mask : ", score);
            printPattern( pattern );
            printf("%.6f before : ", des.getFitness());
            for(int i = 0; i < des.getLength(); i++)
                cout << des.getVal(i);
            cout << endl;
            
            printf("%.6f after  : ", trial.getFitness());
            for(int i = 0; i < trial.getLength(); i++)
                cout << trial.getVal(i);
            cout << endl;
            
            //it = patternList.erase(it);
            //cout << "erasing..." << endl;
            //patternList.push_front( make_pair(score, pattern) );
            //cout << "push_front..." << endl;
            cin.sync();
            cin.get();
            //#endif

            pHash.erase(des.getKey());
            pHash[trial.getKey()] = trial.getFitness();

            des = trial;
            return;
            //if (Chromosome::hit) return;
        }
        if (!evaluated) ++BM_failed;
    }
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
    size_t lastUB = 0;

    for (size_t ub = 1; ub <= mask.size(); ++ub) {

        size_t size = 1;
        Chromosome trial(ell);
        trial = ch;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {
            
            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }
        
        #ifdef DEBUG
        cout << "\n  Try Mask: [";
        for(auto it = mask.begin(); it != next(mask.begin(),size-1); ++it){
            cout << *it << "-";
        }
        cout << "\b]";
        #endif

        //if (isInP(trial)) continue;
        if (isInP(trial)){
            #ifdef DEBUG
            cout << " isInP";
            #endif
            break;
        }
        
        ++RM_failed;
        //2016-10-21
        if (trial.getFitness() > ch.getFitness() - EPSILON) {
        //if (trial.getFitness() >= ch.getFitness()) {
            #ifdef DEBUG
            cout << "\nTaken Mask: [";
            for(auto it = mask.begin(); it != next(mask.begin(),size-1); ++it){
                cout << *it << "-";
            }
            cout << "\b]\n";
            printf(" %.6f before : ", ch.getFitness());
            for(int i = 0; i < ch.getLength(); i++)
                cout << ch.getVal(i);
            cout << endl;
            
            printf(" %.6f after  : ", trial.getFitness());
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

        if (taken) {
            lastUB = ub;
            break;
        }
    }

    if (lastUB != 0) {
        while (mask.size() > lastUB)
            mask.pop_back();
    }

    return taken;

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

        patternList.clear();

        genOrderN();
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
            if (Chromosome::hit) break;
        }

        mergeMasks();

        for (int i=0; i<nCurrent; ++i){
            backMixingO( population[i] );
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
            graph.write(i,j,linkage);
        }
    }


    delete []one;

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
