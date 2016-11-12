/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/

#include <list>
#include <vector>
#include <algorithm>
#include <iterator>
#include <cmath>

#include <iostream>
#include "chromosome.h"
#include "dsmga2.h"
#include "fastcounting.h"
#include "statistics.h"

#include <iomanip>
using namespace std;
#define DEBUG

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

    return termination;

}


bool DSMGA2::foundOptima () {
    return (stFitness.getMax() > population[0].getMaxFitness());
}


void DSMGA2::showStatistics () {

    printf ("Gen:%d  Fitness:(Max/Mean/Min):%f/%f/%f \n ",
            generation, stFitness.getMax (), stFitness.getMean (),
            stFitness.getMin ());
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


void DSMGA2::restrictedMixing(Chromosome& ch) {

    //int r = myRand.uniformInt(0, ell-1);

    //list<int> mask = masks[r];
    
    //2016-11-11 
    list<int> mask;
    int chance = 0; 
    for( pair< list<int>, double > p : sortedMasks ){
        if( chance > ell ) break; 
        mask = p.first;
        //double score = p.second;
        //if (*(mask.begin()) == r){
            //printMaskScore(p);
        //    break;
        //}
    
    
        size_t size = findSize(ch, mask);
        if (size > (size_t)ell/2)
            size = ell/2;
  
        //2016-11-12
        if (mask.size() > size) continue;
        // prune mask to exactly size
        // while (mask.size() > size)
        //    mask.pop_back();

        //2016-11-11 
        //printMaskScore(make_pair(mask, score));

        bool taken = restrictedMixing(ch, mask);
        if (Chromosome::hit) break;

        EQ = true;
        if (taken) {

            genOrderN();

            for (int i=0; i<nCurrent; ++i) {

                if (EQ)
                    backMixingE(ch, mask, population[orderN[i]]);
                else
                    backMixing(ch, mask, population[orderN[i]]);
                if (Chromosome::hit) break;
                
            }
            if (Chromosome::hit) break;
            
            //2016-10-22
            #ifdef DEBUG1
            cout << "population:" << endl;
            for (int i = 0; i < nCurrent; ++i){
                cout << setw(16) << " ";
                for (int j = 0; j < ell; ++j){
                    cout << population[i].getVal(j);
                }
                cout << endl;
            }
            #endif
        }
        break;

    }

}

void DSMGA2::backMixing(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();
        des = trial;
        return;
    }

}

void DSMGA2::backMixingE(Chromosome& source, list<int>& mask, Chromosome& des) {

    Chromosome trial(ell);
    trial = des;
    for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it)
        trial.setVal(*it, source.getVal(*it));

    if (trial.getFitness() > des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        EQ = false;
        des = trial;
        return;
    }

    //2016-10-21
    if (trial.getFitness() > des.getFitness() - EPSILON) {
    //if (trial.getFitness() >= des.getFitness()) {
        pHash.erase(des.getKey());
        pHash[trial.getKey()] = trial.getFitness();

        des = trial;
        return;
    }

}

bool DSMGA2::restrictedMixing(Chromosome& ch, list<int>& mask) {

    bool taken = false;
    size_t lastUB = 0;

    //for (size_t ub = 1; ub <= mask.size(); ++ub) {
    for (size_t ub = mask.size(); ub <= mask.size(); ++ub) {

        size_t size = 1;
        Chromosome trial(ell);
        trial = ch;

        //2016-10-19
        list<int> takenMask;

        for (list<int>::iterator it = mask.begin(); it != mask.end(); ++it) {

            //2016-10-19
            takenMask.push_back(*it);

            trial.flip(*it);

            ++size;
            if (size > ub) break;
        }

        //if (isInP(trial)) continue;
        //if (isInP(trial)) break;
        #ifdef DEBUG
        if (isInP(trial)){
            printMask(takenMask);
            cout << " isInP\n";
        }
        #endif
            //2016-10-22
            #ifdef PRINTPOPULATION 
            cout << "population:" << endl;
            for (int i = 0; i < nCurrent; ++i){
                cout << setw(16) << " ";
                for (int j = 0; j < ell; ++j){
                    cout << population[i].getVal(j);
                }
                cout << endl;
            }
            #endif
            //break;
        //}
        /////////
        #ifdef DEBUG
        cout << " Try Mask:";
        //printMask(takenMask);
        for( const auto &maskPair: sortedMasks ){
           if (maskPair.first == takenMask){
               printMaskScore( maskPair );
               break;
           }
        }

        cout << setw(6) << ch.getFitness() << " before : ";
        for(int i = 0; i < ch.getLength(); i++)
            cout << ch.getVal(i);
        cout << endl;
        
        cout << setw(6) << trial.getFitness() << " after  : ";
        for(int i = 0; i < trial.getLength(); i++)
            cout << trial.getVal(i);
        cout << endl;
        cout << "nfe:" << Chromosome::nfe << ", lsnfe:" << Chromosome::lsnfe << ", hit:" << Chromosome::hit << endl;
        cout << endl ;
        cin.sync();
        cin.get();
        #endif
        ////////
        
        //2016-10-21
        if (trial.getFitness() > ch.getFitness() - EPSILON) {
        //if (trial.getFitness() >= ch.getFitness()) {
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

void DSMGA2::printMask( const list<int> &mask ){
    cout << "[";
    for( const auto &m : mask ){
        cout << m << "-";
    }
    cout << "\b]" << endl;
}

void DSMGA2::printMaskScore( const pair< list<int>, double > &p ){
    list<int> mask = p.first;
    double score = p.second;
    cout << "[";
    for( const auto &m : mask ){
        cout << m << "-";
    }
    cout << "\b] : " << setprecision(10) << score << endl;
}

void DSMGA2::generateRestMask( const list<int> &mask, vector<int> &rest ){
    bool *inMask = new bool[ell]();
    for (auto it = mask.begin(); it != mask.end(); ++it) {
        inMask[ell] = true;
    }
    for (int i = 0; i < ell; ++i) 
        if (!inMask[i])
            rest.push_back(i);
}

double DSMGA2::calculateScore( const list<int> &mask ) {
    if (mask.size() == 1)
        return 0.0;//TODO

    double intraCluster = 0.0;
    for (auto it = mask.begin(); next(it,1) != mask.end(); ++it) {
        for (auto it1 = next(it,1); it1 != mask.end(); ++it1){
            intraCluster += graph( *it, *it1 );
            //cout << "(" << *it << "," <<  *it1 << ")";
        }
    }
    intraCluster /= (mask.size() * (mask.size()-1)/2);

    double interCluster = 0.0;
    vector<int> rest;
    generateRestMask( mask, rest );
    for (auto it = mask.begin(); it != mask.end(); ++it) {
        for (auto it1 = rest.begin(); it1 != rest.end(); ++it1){
            interCluster += graph( *it, *it1 );
            //cout << "(" << *it << "," <<  *it1 << ")";
        }
    }
    interCluster /= (mask.size() * rest.size());


    double score = intraCluster - interCluster;
    //printMaskScore( make_pair( mask, score) );
    return score;
}

bool DSMGA2::maskExist( const list<int> & mask ){
    // Slow version... Need to use map in future
    for( const auto &maskPair: sortedMasks ){
       if (maskPair.first == mask) 
           return true;
    }
    return false;
}

void DSMGA2::sortMasks() {
    for (int i = 0; i < ell; ++i){
        list<int> fullMask = masks[i]; //full length mask, need to split to ILS
        for (int n = 0; n < ell; ++n) {
            list<int> mask = fullMask;
            //printMask(mask);
            mask.sort();
            if (!maskExist( mask )){
                double score = calculateScore(mask);
                sortedMasks.push_back( make_pair( mask, score) );
            }
            else{
                //Replace with higher score ?
            }
            fullMask.pop_back();
        }
    }

    sort( sortedMasks.begin(), sortedMasks.end(),
        [](const pair< list<int>,double > &left, const pair< list<int>,double > &right){
            if (fabs(left.second-right.second)<EPSILON) // float equal comparison
                return left.first.size() > right.first.size();
            else
                return left.second > right.second;
        });

    #ifdef DEBUG 
    for (auto it = sortedMasks.begin(); it != sortedMasks.begin()+50; ++it )
        printMaskScore( *it );
    cin.sync();
    cin.get();
    #endif
}

void DSMGA2::mixing() {

    if (SELECTION)
        selection();

    //* really learn model
    buildFastCounting();
    buildGraph();

    for (int i=0; i<ell; ++i)
        findClique(i, masks[i]);

    //2016-11-11
    sortedMasks.clear();
    sortMasks();

    int repeat = (ell>50)? ell/50: 1;

    for (int k=0; k<repeat; ++k) {

        genOrderN();
        for (int i=0; i<nCurrent; ++i) {
            restrictedMixing(population[orderN[i]]);
            if (Chromosome::hit) break;
        }
        if (Chromosome::hit) break;
    }

    //2016-11-11
    sortedMasks.clear();

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
