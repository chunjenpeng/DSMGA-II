/***************************************************************************
 *   Copyright (C) 2015 Tian-Li Yu and Shih-Huan Hsu                       *
 *   tianliyu@ntu.edu.tw                                                   *
 ***************************************************************************/


#include <math.h>
#include <iostream>
#include <cstdlib>

#include "statistics.h"
#include "dsmga2.h"
#include "global.h"
#include "chromosome.h"

using namespace std;


int
main (int argc, char *argv[]) {


    if (argc != 9) {
        printf ("DSMGA2 ell nInitial function maxGen maxFe repeat display rand_seed\n");
        printf ("function: \n");
        printf ("     ONEMAX:  0\n");
        printf ("     MK    :  1\n");
        printf ("     FTRAP :  2\n");
        printf ("     CYC   :  3\n");
        printf ("     NK    :  4\n");
        printf ("     SPIN  :  5\n");
        printf ("     SAT   :  6\n");

        return -1;
    }

    int ell = atoi (argv[1]); // problem size
    int nInitial = atoi (argv[2]); // initial population size
    int fffff = atoi (argv[3]); // function
    int maxGen = atoi (argv[4]); // max generation
    int maxFe = atoi (argv[5]); // max fe
    int repeat = atoi (argv[6]); // how many time to repeat
    int display = atoi (argv[7]); // display each generation or not
    int rand_seed = atoi (argv[8]);  // rand seed


    if (fffff == 4) {

        char filename[200];
        sprintf(filename, "./NK_Instance/pnk%d_%d_%d_%d", ell, 4, 1, 1);

        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        FILE *fp = fopen(filename, "r");
        loadNKWAProblem(fp, &nkwa);
        fclose(fp);
    }

    if (fffff == 5) {
        char filename[200];
        sprintf(filename, "./SPIN/%d/%d_%d",ell, ell, 1);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSPIN(filename, &mySpinGlassParams);
    }

    if (fffff == 6) {
        char filename[200];
        sprintf(filename, "./SAT/uf%d/uf%d-0%d.cnf", ell, ell, 1);
        if (SHOW_BISECTION) printf("Loading: %s\n", filename);
        loadSAT(filename, &mySAT);
    }


    if (rand_seed != -1)  // time
        myRand.seed((unsigned long)rand_seed);

    int i;

    Statistics stGen, stFE, stLSFE;
    int usedGen;

    int failNum = 0;


    for (i = 0; i < repeat; i++) {

        DSMGA2 ga (ell, nInitial, maxGen, maxFe, fffff);

        if (display == 1)
            usedGen = ga.doIt (true);
        else
            usedGen = ga.doIt (false);


        if (!ga.foundOptima()) {
            failNum++;
            printf ("-");
        } else {
            stFE.record (Chromosome::hitnfe);
            stLSFE.record (Chromosome::lsnfe);
            stGen.record (usedGen);
            printf ("+");
        }

        fflush (NULL);

    }

    cout << endl;
    printf ("\n");
    printf ("%f  %f  %f %d\n", stGen.getMean (), stFE.getMean(), stLSFE.getMean(), failNum);

    if (fffff == 4) freeNKWAProblem(&nkwa);

    return EXIT_SUCCESS;
}
