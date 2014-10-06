/* rc.c - Returns the reverse compliment of a set of sites
 *
 *  MATLAB usage: rcsites = rc(sites)
 *
 * 
 */


#include "mex.h"

#define SITES prhs[0]
#define RCSITES plhs[0]

void mexFunction(   int nlhs,       mxArray *plhs[],
                    int nrhs, const mxArray *prhs[])
{
    // Check for errors
    if (nrhs != 1)
        mexErrMsgTxt("Funciton takes exactly 1 argument1.");
    
    int N = mxGetM(SITES);
    int L = mxGetN(SITES);
    mxChar *sites = (mxChar *) mxGetData(SITES);
    
    int dims[2];
    dims[0] = N;
    dims[1] = L;
    RCSITES = mxCreateCharArray(2, dims);
    mxChar *rcsites = (mxChar *) mxGetData(RCSITES);
    
    char c;
    int i, n;
    for (n=0; n<N; n++) {
        for (i=0; i<L; i++) {
            switch (sites[N*i + n]) {
                case (mxChar) 'A':
                    c = 'T';
                    break;
                case (mxChar) 'C':
                    c = 'G';
                    break;
                case (mxChar) 'G':
                    c = 'C';
                    break;
                case (mxChar) 'T':
                    c = 'A';
                    break;
                case (mxChar) 'a':
                    c = 't';
                    break;
                case (mxChar) 'c':
                    c = 'g';
                    break;
                case (mxChar) 'g':
                    c = 'c';
                    break;
                case (mxChar) 't':
                    c = 'a';
                    break;
                default:
                    c = 'N';
            }
            rcsites[N*(L-i-1) + n] = (mxChar) c;
        }
    }
}