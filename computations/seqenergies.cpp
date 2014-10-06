/* seqenergies.c - computes energy of all valids sites in a sequence
 *
 * MATLAB usage: [energies, positions, orientations] = ...
 *              seqenergies(seq, emat, pos1, pos2)
 *
 *Note: sequences must be passed as doubles with the following correspondences:
 *  1 = A, 2 = C, 3 = G, 4 = T.
 */

#include "mex.h"
#include <math.h>

#define SEQ         prhs[0]
#define EMAT        prhs[1]
#define POS1        prhs[2]
#define POS2        prhs[3]

#define ENERGIES        plhs[0]
#define POSITIONS       plhs[1]
#define ORIENTATIONS    plhs[2]

void mexFunction(   int nlhs,       mxArray *plhs[],
                    int nrhs, const mxArray *prhs[])
{
    // Check for errors
    if (nrhs != 4)
        mexErrMsgTxt("Funciton takes 4 arguments.");
    int position1 = ((int) mxGetScalar(POS1))-1;
    int position2 = ((int) mxGetScalar(POS2))-1;
    if (mxGetM(SEQ) != 1)
        mexErrMsgTxt("Sequence must be a 1xN char matrix.");
    if (mxGetN(EMAT) != 4)
        mexErrMsgTxt("Motif does not have 4 columns.");
    
    // Get seq and site lengths
    int N = mxGetN(SEQ); // Length of sequence
    int L = mxGetM(EMAT); // Length of site
    
    // Get sequence
    double *seq = (double *) mxGetPr(SEQ); // Pointer to sequence
    
    // Get emats
    double *emat = (double *) mxGetPr(EMAT); // Pointer to emat
    double *rcemat = (double *) mxMalloc(4*L*sizeof(double)); //new double[4*L];
    int i, j, k, m, n, x;
    for (i=0; i<L; i++) {
        for (j=0; j<4; j++) {
            m=L-i-1;
            k=4-j-1;
            rcemat[L*k + m] = emat[L*j + i];
        }
    }
    
    // Make sure position of interest is within bounds
    if ((position1 < 0) || (position1 >= L)) {
        mexPrintf("position1 == %d\n", position1);
        mexErrMsgTxt("position1 is out of bounds.");
    }
    if ((position2 < 0) || (position2 >= L)) {
        mexPrintf("position2 == %d\n", position2);
        mexErrMsgTxt("position2 is out of bounds.");
    }
    // Make sure two positions are not the same.
    if (position1==position2)
        mexErrMsgTxt("Two position must be different.");
    
    // Deterimine number of possible sites in the sequence
    int numForwardSites = N-L+1;
    if (numForwardSites < 0)
        numForwardSites = 0;
    
    // Create vector of energies to be returned 
    ENERGIES = mxCreateDoubleMatrix(2*numForwardSites,1,mxREAL);
    double *energies = (double *) mxGetPr(ENERGIES); 
    
    // Create vector of positions to be returned
    POSITIONS = mxCreateDoubleMatrix(2*numForwardSites,1,mxREAL);
    double *positions = (double *) mxGetPr(POSITIONS);
    
    // Create vector of orientations to be returned
    // Orientations also contains info on which sites
    // have which base at the positionOfInterest
    ORIENTATIONS = mxCreateDoubleMatrix(2*numForwardSites,1,mxREAL);
    double *orientations = (double *) mxGetPr(ORIENTATIONS);
    
    // Initialize arrays, assuming all sites are valid
    for (n=0; n<numForwardSites; n++) {
        orientations[n] = 1;
        orientations[numForwardSites+n] = -1;
        positions[n] = n+1;
        positions[numForwardSites+n] = n+1;
        energies[n] = 0.0;
        energies[numForwardSites+n] = 0.0;
    }
    
    // Compute which sites are invalid
    int badstart, badstop;
    bool valid = true;
    int nfwd, nrc;
    for (i=0; i<N; i++) {
        if (seq[i] == 0.)    {
            badstart = i-L+1;
            if (badstart < 0)
                badstart = 0;
            badstop = i+1;
            if (badstop > numForwardSites)
                badstop = numForwardSites;
            for (n=badstart;n<badstop;n++) {
                orientations[n] = 0;
                orientations[numForwardSites+n] = 0;
            }
        }
    }
    
    // Compute energies of all valid sites
    double energy, rcenergy;
    bool match, rcmatch;
    int* bases, rcbases;
    for (n=0;n<numForwardSites;n++) {
        if (orientations[n] != 0) { // Means orientaions[numForwardSites+n] != 0 too.
            energy = 0.0;
            rcenergy = 0.0;
            
            // Compute energy matrix contributions
            for (i=0;i<L;i++)   {
                x = L*(((int)seq[n+i])-1)+i;
                energy += emat[x];
                rcenergy += rcemat[x];
            }
            
            // Record energies
            energies[n] = energy;
            energies[numForwardSites+n] = rcenergy;
            
            // Record bases at given site
            orientations[n] = 
                4*(((int)seq[n+position1])-1)+(int)seq[n+position2];
            orientations[numForwardSites+n] =
                -1*(4*(4-((int)seq[n+(L-(position1+1))])) + (5-((int)seq[n+(L-(position2+1))])));
        }
    }  
        
    // Free memory. Matlab says it does this, but keeps crashing if I don't.
    //delete[] rcemat;
    mxFree(rcemat);
}


