/*
 * Generating InvQ from the conditioned distribution 
 */


#include <R.h>
#include <Rmath.h>

#include "myrwishart.h"
#include "matrmult.h"
#include "cholesky.h"

#include "structures.h"

void gibbs_InvQ(struct str_state* last,  // OUT+IN last known values of generated parameters
                struct str_param* param, // IN hyperparameters
                int* spec,            // [7]                class-specific parameters
                                        // order:   [0] gamma, 
                                        //          [1] tau, 
                                        //          [2] beta, 
                                        //          [3] rantau, 
                                        //          [4] mu, 
                                        //          [5] InvQ, 
                                        //          [6] InvSigma
                int* calc,            // [9]                calculable parameters
                                        // order:   [0]   pUik_nonb, 
                                        //          [1-3] Q, detQ, detInvQ, 
                                        //          [4-6] Sigma, detSigma, detInvSigma, 
                                        //          [7,8] sdSigma, corSigma
                int* dims,            // [23]       the length of subarray that corresponds to one state (disected by various parameters)
                                        // order  [0-2]   w, pUik, U,
                                        //        [3-6]   gamma, min, max, latent,
                                        //        [7-8]   beta, tau,
                                        //        [9-10]  mu, rantau,
                                        //        [11-12] InvSigma, InvQ,
                                        //        [13]    b,
                                        //        [14]    pUik_nonb,
                                        //        [15-16] Sigma, Q,
                                        //        [17-20] det(Sigma,InvSigma,Q,InvQ),
                                        //        [21-22] sdSigma, corSigma,
                int* dimswithK,       // [23]       the length of subarray that corresponds to one state 
                                        //            (disected by various parameters, also multiplication by K incorporated when such parameters is class-specific)
                                        // order  [0-2]   w, pUik, U,
                                        //        [3-6]   gamma, min, max, latent,
                                        //        [7-8]   beta, tau,
                                        //        [9-10]  mu, rantau,
                                        //        [11-12] InvSigma, InvQ,
                                        //        [13]    b,
                                        //        [14]    pUik_nonb,
                                        //        [15-16] Sigma, Q,
                                        //        [17-20] det(Sigma,InvSigma,Q,InvQ),
                                        //        [21-22] sdSigma, corSigma,
                                        
                int* K,                 // IN  [1] number of classes
                int* totnran,           // IN [1] total number of random regressors
                int* i,                 // IN [1] iteration number
                double* Q,              // [(B+M)*totnran*(totnran+1)/2(* K)] Prior scale matrix of InvSigma
                double* detQ,           // [(B+M)(* K)] determinant of Q
                double* detInvQ         // [(B+M)(* K)] determinant of InvQ
                )    
{ 
  /** Declarations **/
  int j, k;                 // looping indeces
  int coord;                // coordinate for matrices as arrays
  double Smatrix[dims[12]]; // Scale matrix of Wishart distribution
  double M[dims[12]];       // sum of InvSigma and param$V
  double chol[dims[12]];    // cholesky decomposition of matrix M
  double df;                // degrees of freedom for Wishart distribution
  
  // Is InvQ class specific?
  if(spec[5]){
    // and so is InvSigma...
    
    for(k = 0; k < *K; k++){
      // find the scale matrix
      for(j = 0; j < dims[12]; j++){
        M[j] = (*last).InvSigma[j + k*dims[12]] + (*param).InvV[j];
      }

      justcholesky(M, chol, totnran);
      choltoinv(chol, Smatrix, totnran);
      
      // cholesky decomposition of Smatrix for generating from wishart
      justcholesky(Smatrix, chol, totnran);
        
      df = *((*param).nu_0) + *((*param).nu_1);
      
      my_rwishart((*last).InvQ + k*dims[12], chol, &df, totnran);
      
      // calculate possible parameters
      if( calc[1] || calc[2] || calc[3]){
        justcholesky((*last).InvQ + k*dims[12], chol, totnran);
      }
      
      if(calc[1]){
        choltoinv(chol, Q + *i * dimswithK[12] + k*dims[12], totnran);
      }
      
      if(calc[2]){
        coord = *i * *K + k;
        detchol(chol, detQ + coord, totnran);
        detQ[coord] = 1/detQ[coord];
      }
      
      if(calc[3]){
        coord = *i * *K + k;
        detchol(chol, detInvQ + coord, totnran);
      }
      
    }
  }else{
    // InvQ is just one for all classes
    
    // Is InvSigma class-specific?
    if(spec[6]){
      // inicialization of Smatrix
      for(j = 0; j < dims[12]; j++){
        M[j] = (*param).InvV[j];
      }
      // add InvSigma for each class
      for(k = 0; k < *K; k++){
        for(j = 0; j < dims[12]; j++){
          M[j] += (*last).InvSigma[j + k * dims[12]];
        }
      }
      
      df = *K * *((*param).nu_0) + *((*param).nu_1);
    }else{
      
      for(j = 0; j < dims[12]; j++){
        M[j] = (*param).InvV[j] + (*last).InvSigma[j];
      }
      
      df = *((*param).nu_0) + *((*param).nu_1);
    } // end of if(spec["InvSigma"])
    
    // finding the inverse
    justcholesky(M, chol, totnran);
    choltoinv(chol, Smatrix, totnran);
    // cholesky decomposition of Smatrix for generating from wishart
    justcholesky(Smatrix, chol, totnran);
    
    my_rwishart((*last).InvQ, chol, &df, totnran);
      
    
    // calculate possible parameters
    if( calc[1] || calc[2] || calc[3]){
      justcholesky((*last).InvQ, chol, totnran);
    }
    
    if(calc[1]){
      choltoinv(chol, Q + *i * dimswithK[12], totnran);
    }
    
    if(calc[2]){
      detchol(chol, detQ + *i, totnran);
      detQ[*i] = 1/detQ[*i];
    }
    
    if(calc[3]){
      detchol(chol, detInvQ + *i, totnran);
    }
  } // end of if(spec["InvQ"])
  
} // end of void 
