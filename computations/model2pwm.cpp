/* model2pwm.c - computes pwm given emat and cutoff
 * If model has multiple cutoffs, only first is used
 *
 * MATLAB usage: [pwm, energies] = model2pwm(model, numTrials);
 *
 */

#include "mex.h"
#include <cstdlib>
#include <ctime>
#include <iostream>

#define MODEL       prhs[0]
#define NUMTRIALS   prhs[1]
#define PWM         plhs[0]
#define ENERGIES    plhs[1]

void mexFunction(   int nlhs,       mxArray *plhs[],
                    int nrhs, const mxArray *prhs[])
{
    // Check for errors
    int N;
    if (nrhs == 2)
        N = (int) mxGetScalar(NUMTRIALS); 
    else if (nrhs == 1)
        N = 1000;
    else
        mexErrMsgTxt("Funciton takes 1 or 2 arguments."); 
    
    // Get model.cutoff(1)
    mxArray *CUTOFF = mxGetField(MODEL, 0, "cutoff");
    if (!CUTOFF)
        mexErrMsgTxt("Input 'model' does not have field 'cutoff'");
    double cutoff = mxGetScalar(CUTOFF);

    // Get model.emat
    mxArray *EMAT = mxGetField(MODEL, 0, "emat");
    if (!EMAT)
        mexErrMsgTxt("Input 'model' does not have field 'emat'");
    if (mxGetN(EMAT) != 4)
        mexErrMsgTxt("Emat does not have 4 columns.");
    double *emat = (double *) mxGetPr(EMAT);
    unsigned L = mxGetM(EMAT); 

    // Create pwm
    PWM = mxCreateDoubleMatrix(L,4,mxREAL);
    double *pwm = (double *) mxGetPr(PWM);
    
    // Create emat
    ENERGIES = mxCreateDoubleMatrix(N,1,mxREAL);
    double *energies = (double *) mxGetPr(ENERGIES); 
    
    // Declare variables
    unsigned i,n,b,bold,bnew,bmin;
    double energy, energyold, x, xmin;
    unsigned *site = (unsigned*) mxMalloc(L*sizeof(unsigned)); //new unsigned[L];
    unsigned *siteold = (unsigned*) mxMalloc(L*sizeof(unsigned)); //new unsigned[L];
    
    // Initialize pwm
    for (i=0;i<4*L;i++) {
        pwm[i] = 0;
    }
    
    // Find consensus site
    energy = 0.0;
    for (i=0;i<L;i++) {
        bmin = 0;
        xmin = emat[L*bmin + i];
        for (b=1;b<4;b++) {
            x = emat[L*b + i];
            if (x < xmin)   {
                bmin = b;
                xmin = x;
            }   
        }
        site[i] = bmin;
        energy += xmin;
    }
    
    // Make sure energy of consensus site is <= than cutoff
    if (energy >= cutoff) {
        cutoff = energy;
        //mexErrMsgTxt("ERROR: Consensus site energy is greater or equal to cutoff.");
    }
    // Do MCMC
    srand((unsigned)time(0));
    for (n=0;n<N;n++) {
        // Record old site and energy
        for (i=0;i<L;i++)
            siteold[i] = site[i];
        energyold = energy;
        
        // Choose random position
        i = rand()%L;
        bold = siteold[i];
        
        // Choose random new base
        bnew = bold;
        while (bnew == bold)
            bnew = rand()%4;
        site[i] = bnew;
        
        // Compute new energy
        energy = energyold - emat[L*bold+i] + emat[L*bnew+i];
        
        // Reject step if energy is above cutoff
        if (energy >= cutoff) {
            site[i] = bold;
            energy = energyold;
        }
            
        // Record energy and increment pwm
        energies[n] = energy;
        for (i=0;i<L;i++)
            pwm[L*site[i] + i] += 1.0/double(N);
    }
    
    mxFree(site); //delete[] site;
    mxFree(siteold); //delete[] siteold;
}
