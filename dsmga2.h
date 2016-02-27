/*
 * dsmga2.h
 *
 *  Created on: May 2, 2011
 *      Author: tianliyu
 */

#ifndef _DSMGA2_H_
#define _DSMGA2_H_

#include <list>
#include "chromosome.h"
#include "statistics.h"
#include "trimatrix.h"
#include "doublelinkedlistarray.h"
#include "fastcounting.h"

class DSMGA2 {
public:
    DSMGA2 (int n_ell, int n_nInitial, int n_maxGen, int n_maxFe, int fffff);

    ~DSMGA2 ();

    void selection ();
    /** tournament selection without replacement*/
    void tournamentSelection();

    void oneRun (bool output = true);
    int doIt (bool output = true);

    void buildGraph ();
    void mixing ();
    void restrictedMixing(Chromosome&);
    bool restrictedMixing(Chromosome& ch, list<int>& mask);
    void backMixing(Chromosome& source, list<int>& mask, Chromosome& des);
    void backMixingE(Chromosome& source, list<int>& mask, Chromosome& des);

    bool shouldTerminate ();

    bool foundOptima ();

    int getGeneration () const {
        return generation;
    }

    bool isInP(const Chromosome& ) const;
    void genOrderN();
    void genOrderELL();

    void showStatistics ();

    bool isSteadyState ();

//protected:
public:

    int ell;                                  // chromosome length
    int nCurrent;                             // population size
    bool EQ;
    unordered_map<unsigned long, double> pHash; // to check if a chromosome is in the population


    list<int> *masks;
    int *selectionIndex;
    int *orderN;                             // for random order
    int *orderELL;                             // for random order
    int selectionPressure;
    int maxGen;
    int maxFe;
    int repeat;
    int generation;
    int bestIndex;

    Chromosome* population;
    FastCounting* fastCounting;

    TriMatrix<double> graph;

    double previousFitnessMean;
    Statistics stFitness;

    // methods
    double computeMI(double, double, double, double) const;


    void findClique(int startNode, list<int>& result);

    void buildFastCounting();
    int countXOR(int, int) const;
    int countOne(int) const;

    size_t findSize(Chromosome&, list<int>&) const;
    size_t findSize(Chromosome&, list<int>&, Chromosome&) const;


};


#endif /* _DSMGA2_H_ */
