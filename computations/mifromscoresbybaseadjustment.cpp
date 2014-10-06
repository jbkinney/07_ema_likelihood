#include "mex.h"
#include "mi.h"

#define ADJUSTMENTS         prhs[0]
#define SCORESBYBASE        prhs[1]
#define CUTOFF              prhs[2]
#define REGIONTYPES         prhs[3]
#define NUMTYPES            prhs[4]

#define MI                  plhs[0]
#define CZX                 plhs[1]


// Main function
void mexFunction(   int nlhs,       mxArray *plhs[],
                    int nrhs, const mxArray *prhs[])
{
    int i, j, m, n, r;
    double inc;
    double* czx;
    double score, minscore;
    
    double *adjustments = (double *) mxGetPr(ADJUSTMENTS);
    if (mxGetNumberOfElements(ADJUSTMENTS) != 16)
        mexErrMsgTxt("ADJUSTMENTS does not have 16 elements");
    double *scoresByBase = (double *) mxGetPr(SCORESBYBASE);
    if (mxGetN(SCORESBYBASE) != 16)
        mexErrMsgTxt("SCORESBYBASE does not have 16 cols");
    int numRegions = mxGetM(SCORESBYBASE);
    double cutoff = mxGetScalar(CUTOFF);
    double *regionTypes = (double *) mxGetPr(REGIONTYPES);
    if (mxGetNumberOfElements(REGIONTYPES) != numRegions)
        mexErrMsgTxt("REGIONTYPES does not have same number of elements as SCORESBYBASE has rows");
    int numTypes = (int) mxGetScalar(NUMTYPES);
    
    
    MI = mxCreateDoubleMatrix(1,1,mxREAL);
    double* mi = (double*) mxGetPr(MI);
    
    // Initialize joint and marginal distributions
    CZX = mxCreateDoubleMatrix(numTypes,2,mxREAL);
    czx = (double *) mxGetPr(CZX);
    for (i=0;i<2*numTypes;i++) {
        czx[i] = 0.0;
    }
    
    // Compute joint distribution given adjustments
    for (n=0;n<numRegions;n++) {
        minscore = scoresByBase[0*numRegions+n] + adjustments[0];
        for (i=1;i<16;i++) {
            score =  scoresByBase[i*numRegions+n] + adjustments[i];
            if (score < minscore)
                minscore = score;
        }
        m = (minscore < cutoff) ? 1 : 0;
        r = (int) regionTypes[n];
        czx[r+numTypes*m] += 1.0;
    }
    
    mi[0] = normalLogLikelihood(czx, numTypes, 2);
}