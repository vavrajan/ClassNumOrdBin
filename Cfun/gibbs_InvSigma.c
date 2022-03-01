/*
 * Generating InvSigma from the conditioned distribution 
 */


#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>

#include "myrwishart.h"
#include "matrmult.h"
#include "cholesky.h"

#include "structures.h"

void gibbs_InvSigma(struct str_state* last,  // OUT+IN last known values of generated parameters
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
                                            
                    int* K,                 // IN [1] number of classes
                    int* n,                 // IN [1] number of subjects
                    int* totnran,           // IN [1] total number of random regressors
                    int* i,                 // IN [1] iteration number
                    int* nUk,               // IN [K] number of subjects     in classes 1, ..., K
                    double* Sigma,          // [(B+M)*totnran*(totnran+1)/2(* K)] Prior variance matrix of random effects
                    double* detSigma,       // [(B+M)(* K)] determinant of Sigma
                    double* detInvSigma,    // [(B+M)(* K)] determinant of InvSigma
                    double* sdSigma,        // [(B+M)*totnran(* K)] Prior standard deviations of random effects
                    double* corSigma       // [(B+M)*totnran*(totnran-1)/2(* K)] Prior correlation matrix of random effects
                    )    
{ 
  /** Declarations **/
  int j, k, l, m;           // looping indeces
  int coord;                // coordinate for matrices as arrays
  double Smatrix[dims[11]]; // Scale matrix of Wishart distribution
  double M[dims[11]];       // sum of InvSigma and param$V
  double chol[dims[11]];    // cholesky decomposition of matrix M
  double df;                // degrees of freedom for Wishart distribution
  double* pInvQ;            // pointer to InvQ
  double* pmu;              // pointer to mu
  double Sig[dims[11]];     // auxiliary matrix Sigma
  double* pSigma;           // pointer to Sigma  
  
  /*
  for(j = 0; j < dimswithK[11]; j++){
    printf("%d = %f\n", j, (*last).InvSigma[j]);
  }
  */
  
  // Is InvSigma class specific?
  if(spec[6]){
    
    for(k = 0; k < *K; k++){
      /* Setting pointers */
      // InvQ
      pInvQ = (*last).InvQ;
      if(spec[5]){
        pInvQ += k*dims[12];
      }
      // mu
      pmu = (*last).mu;
      if(spec[4]){
        pmu += k*dims[9];
        }
      
      // find the scale matrix
      for(l = 0; l < *totnran; l++){
        for(j = l; j < *totnran; j++){
          coord = l + j*(j+1)/2;
          // inicialization
          M[coord] = pInvQ[coord];
          // add (b-mu)T(b-mu)
          for(m = 0; m < *n; m++){
            // add only if this subject lies in category k
            if((*last).U[m] == k){
              M[coord] += ((*last).b[m + *n * l] - pmu[l]) * ((*last).b[m + *n * j] - pmu[j]);
            }
          }
        }
      }
     
      // positive definite matrix M now needs to be inversed
      justcholesky(M, chol, totnran);
      choltoinv(chol, Smatrix, totnran);
      
      // cholesky decomposition of Smatrix for generating from wishart
      justcholesky(Smatrix, chol, totnran);
      df = nUk[k] + *((*param).nu_0);
      
      /*
        printf("%p\n", (*last).InvSigma + k*dims[11]);
        printf("%p\n", &((*last).InvSigma[k*dims[11]]));
      */
      
      // generating from Wishart distribution
      my_rwishart((*last).InvSigma + k*dims[11], chol, &df, totnran);

      // calculate possible parameters
      // we will surely need cholesky decomposition
      justcholesky((*last).InvSigma + k*dims[11], chol, totnran);
      detchol(chol, (*last).detInvSigma + k, totnran);
      
      
      if(calc[4]){
        choltoinv(chol, Sigma + *i * dimswithK[15] + k*dims[15], totnran);
      }
      
      if(calc[5]){
        detSigma[*i * *K + k] = 1/(*last).detInvSigma[k];
      }
      
      if(calc[6]){
        detInvSigma[*i * *K + k] = (*last).detInvSigma[k];
      }
      
      if(calc[7] || calc[8]){
        // setting pointer to Sigma matrix
        if(calc[4]){
          pSigma = Sigma + *i * dimswithK[15] + k*dims[15];
        }else{ // inverse was not calculated yet
          choltoinv(chol, Sig, totnran);
          pSigma = Sig;
          
        }
        
        /*
        for(l = 0; l < *totnran; l++){
          for(j = l; j < *totnran; j++){
            printf("%f, ", Sig[l + j*(j+1)/2]);
          }
          printf("\n");
        }
        */
        
        if(calc[7]){
          for(j = 0; j < *totnran; j++){
            sdSigma[*i * dimswithK[21] + k*dims[21] + j] = sqrt(pSigma[j*(j+3)/2]);
          }
        }
        
        
        
        if(calc[8]){
          for(l = 0; l < *totnran; l++){
            for(j = l+1; j < *totnran; j++){
              corSigma[*i * dimswithK[22] + k*dims[22] + l + j*(j-1)/2] = pSigma[l + j*(j+1)/2] / sqrt(pSigma[j*(j+3)/2] * pSigma[l*(l+3)/2]); 
            }
          }
        }
      }
    } // end of for k
    
    //---------------------------------------------------------
  }else{
    // InvSigma is just one for all classes --> so must be InvQ
    pInvQ = (*last).InvQ;
    
    // find the scale matrix
    for(l = 0; l < *totnran; l++){
      for(j = l; j < *totnran; j++){
        coord = l + j*(j+1)/2;
        // inicialization
        M[coord] = pInvQ[coord];
        // add (b-mu)T(b-mu)
        for(m = 0; m < *n; m++){
          // add always but mu may depend on this subjects class k
          pmu = (*last).mu;
          if(spec[4]){
            pmu += (*last).U[m] * dims[9];
          }
          
          M[coord] += ((*last).b[m + *n * l] - pmu[l])* ((*last).b[m + *n * j] - pmu[j]);
        }
      }
    }
    
    
    // positive definite matrix M now needs to be inversed
    justcholesky(M, chol, totnran);
    choltoinv(chol, Smatrix, totnran);
    
    // cholesky decomposition of Smatrix for generating from wishart
    justcholesky(Smatrix, chol, totnran);
    df = *n + *((*param).nu_0);
    
    // generating from Wishart distribution
    my_rwishart((*last).InvSigma, chol, &df, totnran);
    
    // calculate possible parameters
    // we will surely need cholesky decomposition
    justcholesky((*last).InvSigma, chol, totnran);
    detchol(chol, (*last).detInvSigma, totnran);
    
    if(calc[4]){
      choltoinv(chol, Sigma + *i * dimswithK[15], totnran);
    }
    
    if(calc[5]){
      detSigma[*i] = 1/(*last).detInvSigma[0]; // [0] needs to be there to not to return adress but the value pointed by this adress
    }
    
    if(calc[6]){
      detInvSigma[*i] = (*last).detInvSigma[0]; // [0] needs to be there to not to return adress but the value pointed by this adress
    }
    
    if(calc[7] || calc[8]){
      // setting pointer to Sigma matrix
      if(calc[4]){
        pSigma = Sigma + *i * dimswithK[15];
      }else{ // inverse was not calculated yet
        choltoinv(chol, Sig, totnran);
        pSigma = Sig;
      }
      
      if(calc[7]){
        for(j = 0; j < *totnran; j++){
          sdSigma[*i * dimswithK[21] + j] = sqrt(pSigma[j*(j+3)/2]);
        }
      }
      
      if(calc[8]){
        for(l = 0; l < *totnran; l++){
          for(j = l+1; j < *totnran; j++){
            corSigma[*i * dimswithK[22] + l + j*(j-1)/2] = pSigma[l + j*(j+1)/2] / sqrt(pSigma[l*(l+3)/2] * pSigma[j*(j+3)/2]); 
          }
        }
      }
    }
  } // end of if(spec["InvSigma"])
  
} // end of void 
