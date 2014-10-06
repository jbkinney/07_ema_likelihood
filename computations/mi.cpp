#include <math.h>
#include "mex.h"
#include "mi.h"

////////////////////////////////////////////////////
// Functions for computing normalized log likelihood
////////////////////////////////////////////////////
// Log gamma funciton
double gammaln(const double x) {
    double f,z;
    
    if (x >= .99999) {
        f = .5*log(2*M_PI) + (x-.5)*log(x) - x + 
            (1./12.)*pow(x,-1.) - (1./360.)*pow(x,-3.) +
            (1./1260.)*pow(x,-5.) - (1./1680.)*pow(x,-7.);
    }
    else    {
        mexErrMsgTxt("Error: Trying to compute gammaln(x) for x < 1");
    }
    return f;    
}

// Compute raw log likelihood per element
double perDatumLogLikelihood(const double *czx, int numz, int numx) {
    double cx;
    int z,x;
    double f,N;
    
    N = 0.;
    f = 0.;
    for (x=0;x<numx;x++)    {
        cx = 0.; 
        for (z=0;z<numz;z++) {
            cx += czx[z+x*numz];
            f += gammaln(czx[z+x*numz]+1);
        }
        N += cx;
        f -= gammaln(cx + numz);
    }
    return f/N;
}

// Compute normalized log likelihood
double normalLogLikelihood(const double *czx, int numz, int numx) { 
    mxArray *CZ = mxCreateDoubleMatrix(numz,1,mxREAL);
    double *cz = (double *) mxGetPr(CZ);
    //cz = new double[numz];

    mxArray *CX = mxCreateDoubleMatrix(numx,1,mxREAL);
    double *cx = (double *) mxGetPr(CX);
    //cx = new double[numx];
    
    mxArray *DZX = mxCreateDoubleMatrix(numz,numx,mxREAL);
    double *dzx = (double *) mxGetPr(DZX);
    //dzx = new double[numz*numx];
    
    int x,z;
    double N = 0.;
    
    for (x=0;x<numx;x++) {
        cx[x] = 0.;
    }
    for (z=0;z<numz;z++) {
        cz[z] = 0.;
    }
    for (x=0;x<numx;x++) {
        for (z=0;z<numz;z++) {
            cz[z] += czx[z+x*numz];
            cx[x] += czx[z+x*numz];
            N += czx[z+x*numz];
        }
    }
    for (x=0;x<numx;x++) {
        for (z=0;z<numz;z++) {
            dzx[z+x*numz] = cz[z]*cx[x]/N;
        }
    }
    
    return perDatumLogLikelihood(czx,numz,numx) - perDatumLogLikelihood(dzx,numz,numx);
    
    //delete[] cz;
    //delete[] cx;
    //delete[] dzx;
}


////////////////////////////////////////////////////
// Functions for computing entropy and mutual info
////////////////////////////////////////////////////

// Useful in entropy and mutual info computations
double xlnx(const double p) {
    double f;
    if (p > 0.0)
        f = p*log(p);
    else 
        f = 0.0;
    return f;
}


// Function to compute (naive) entropy
double entropy(const double* prob, const int probLength) {
    int i;
    double H = 0.0;
    for (i=0;i<probLength;i++)
        H -= xlnx(prob[i]);
    return H;       
}



// Function to compute (naive) mutual information
double mutualinfo(const double *cij, const int numi, const int numj) {
    int i,j,k;
    double I;
    
    // Compute marginal distributions
    double *pi = (double *) mxMalloc(numi*sizeof(double)); //new double[numi];
    double *pj = (double *) mxMalloc(numj*sizeof(double)); //new double[numj];
    double *pij = (double *) mxMalloc(numi*numj*sizeof(double)); //new double[numi*numj];
    double N =0.;
    for (i=0; i<numi; i++)  {
        pi[i] = 0.0;
    }
    for (j=0; j<numj; j++)  {
        pj[j] = 0.0;
    }
    for (k=0; k<numi*numj; k++) {
        N += cij[k];
        pij[k] = 0.;
    }
    
    for (i=0; i<numi; i++) {
        for (j=0; j<numj; j++) {
            pij[i+j*numi] = cij[i+j*numi]/N;
            pi[i] += pij[i+j*numi];
            pj[j] += pij[i+j*numi];
        }
    }
    
    // Compute MI using entropy of joint and marginal distributions
    I = entropy(pi,numi) + entropy(pj,numj) - entropy(pij, numi*numj);
    mxFree(pi); //delete[] pi;
    mxFree(pj); //delete[] pj;
    mxFree(pij); //delete[] pij;
    return I;
}

////////////////////////////////////////////////////
// Other stuff
////////////////////////////////////////////////////

/*
// Log factorial funciton
double lf(const double x) {
    double z,f;
    
    // Use Stirling's approx if x > 10
    if (x > 0) {
        z = x + 1.0;
        f = .5*log(2*M_PI) + (z-.5)*log(z) - z + 
        (1./12.)*pow(z,-1.) - (1./360.)*pow(z,-3.) +
        (1./1260.)*pow(z,-5.) - (1./1680.)*pow(z,-7.);
    }
    else
        f = 0.0;
    return f;
}



// Function to compute correction to mutual info using uniform 
// prior on space of error models
// j indexes over predictions, i over data;
double adjustedmutualinfo(const double *cij, const int numi, const int numj) {
    int i,j,k;
    double I;
    double N = 0.0;
    
    // Compute marginal distributions
    double *pi = new double[numi];
    double *pj = new double[numj];
    double *pij = new double[numi*numj];
    double *cj = new double[numj];
    for (k=0; k<numi*numj; k++) {
        N += cij[k];
    }
    for (i=0; i<numi; i++)  {
        pi[i] = 0.0;
    }
    for (j=0; j<numj; j++) { 
        pj[j] = 0.0;
        cj[j] = 0.0;
    }
    for (i=0; i<numi; i++) {
        for (j=0; j<numj; j++) {
            pij[i+j*numi] = cij[i+j*numi]/N;
            pi[i] += cij[i+j*numi]/N;
            pj[j] += cij[i+j*numi]/N;
            cj[j] += cij[i+j*numi];  
        }
    }
    
    // Compute MI using entropy of joint and marginal distributions
    I = entropy(pi,numi) + entropy(pj,numj) - entropy(pij, numi*numj);
    
    double A = 0.;
    double B = 0.;
    double C = 0.;
    double D = 0.;
    double E = 0.;
    double Delta = 0.;
    A = -1.*double(numj)*lf(double(numi-1))/N; 
    for(j=0;j<numj;j++) {
        B += lf(numi-1+cj[j])/N;
    }
    for(j=0;j<numj;j++) {
        C -= xlnx(cj[j])/N;
    }
    for(k=0; k<numi*numj; k++) {
        D -= lf(cij[k])/N;
    }
    for(k=0; k<numi*numj; k++)  {  
        E += xlnx(cij[k])/N;
    }
    //mexPrintf("I == %f\n", I);
    //mexPrintf("A == %f\n", A);
    //mexPrintf("B == %f\n", B);
    //mexPrintf("C == %f\n", C);
    //mexPrintf("D == %f\n", D);
    //mexPrintf("E == %f\n", E);
    
    delete[] pi;
    delete[] pj;
    delete[] pij;
    delete[] cj;
    
    Delta = A+B+C+D+E;
    return I-Delta;
}
 **/