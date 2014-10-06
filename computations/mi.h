#include <math.h>
#include "mex.h"

// Log gamma funciton
double gammaln(const double x);

// Compute raw log likelihood per element
double perDatumLogLikelihood(const double *czx, int numz, int numx);

// Compute normalized log likelihood
double normalLogLikelihood(const double *czx, int numz, int numx);

// Useful in entropy and mutual info computations
double xlnx(const double p);

// Function to compute (naive) entropy
double entropy(const double* prob, const int probLength);

// Function to compute (naive) mutual information
double mutualinfo(const double *cij, const int numi, const int numj);

