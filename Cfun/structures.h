#ifndef STRUCTURES_H
#define STRUCTURES_H

struct str_param {
  double* nu_0;     // [1]    prior df for InvSigma
  double* nu_1;     // [1]    prior df for InvQ  
  double* delta;    // [1]    prior parameter of Dirichlet distribution
  double* a;        // [1]    prior shape of tau
  double* b;        // [1]    prior rate  of tau
  double* aran;     // [1]    prior shape of rantau
  double* bran;     // [1]    prior rate  of rantau
  double* InvV;     // [totnran*(totnran+1)/2] prior inverse scale matrix for InvQ
  double* ranmu0;   // [totnran] prior mean of mu
  double* gamma1;   // [nY[1 ... Ords]] value of first gamma parameter (for every Ord variable)
  double* gammaBin; // [1]    value of gamma parameter for Bin variables (all the same)
  double* fixD;     // [totnfix] prior variance of fixed beta parameters
  double* fixmu;    // [totnfix] prior mean value of fixed beta parameters
};

struct str_state {
  double* w;        // [K]                    Class probabilities
  int*    U;        // [n]                    To which class does i-th subject belong to
  double* pUik;     // [n * K]                Probability of i-th subject belong to each class
  double* latent;   // [(nY[1]+nY[2])*N]      Latent variables
  double* gamma;    // [sum(ncat-2)(* K)]     Division points 
  double* beta;     // [totnfix(* K)]         Beta parameters for fixed effects
  double* b;        // [totnran * n]          Random effects
  double* tau;      // [nY[0](* K)]           Precision parameter of numerical variables
  double* mu;       // [totnran(* K)]         Mean values of random effects
  double* rantau;   // [totnran(* K)]         Prior precision parameters of mu
  double* InvSigma; // [totnran*(totnran+1)/2(* K)] Precision matrix of random effects
  double* InvQ;     // [totnran*(totnran+1)/2(* K)] Prior inverse scale matrix of InvSigma
  double* detInvSigma; // [1 (*K)]            determinant of InvSigma (calculable parameter)
};

#endif
