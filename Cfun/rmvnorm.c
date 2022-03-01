/*
 * Generating from multivariate normal distribution of dimension p
 * using mean vector mu and INVERSE variance matrix (Precision matrix) InvSigma
 */


#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "matrmult.h"
#include "cholesky.h"

void rmvnorm(double* out,         // [p] generated vector
             int*    p,           // [1] dimension
             double* mu,          // [p] mean value
             double* InvSigma)    // [p*(p+1)/2] symmetric precision matrix
{
  double chol [*p * (*p+1)/2];   // Cholesky decomposition
  double gen [*p];              // N(0,1) values
  int ii;                      // looping index
  
  // Generate X ~ N(mu, Sigma) as follows
  // X = C^{-1} * Y + mu
  // InvSigma = C^{T} * C -=- Cholesky decomposition (chol upper triangular)
  
  // first do Cholesky decomposition
  justcholesky(InvSigma, chol, p);
  
  // generate Y ~ N(0,I)
  for(ii = 0; ii < *p; ii++){
    gen[ii] = rnorm(0.0, 1.0);
  }
  
  // now backsolve 
  backsolve2(chol, gen, p, out);
  
  // add mu
  for(ii = 0; ii < *p; ii++){
    out[ii] += mu[ii];
  }
  
}
