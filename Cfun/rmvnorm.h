#ifndef RMVNORM_H
#define RMVNORM_H

void rmvnorm(double* out,         // [p] generated vector
             int*    p,           // [1] dimension
             double* mu,          // [p] mean value
             double* InvSigma);   // [p*(p+1)/2] symmetric precision matrix
 
#endif