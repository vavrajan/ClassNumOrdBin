/*
 * Generating tau from the conditioned distribution 
 */


#include <R.h>
#include <Rmath.h>

#include "structures.h"

void gibbs_tau( struct str_state* last,  // OUT+IN last known values of generated parameters
                struct str_param* param, // IN hyperparameters
                double* Y,            // IN [*N * totnY]
                int* spec,            // [7]                class-specific parameters
                                        // order:   [0] gamma, 
                                        //          [1] tau, 
                                        //          [2] beta, 
                                        //          [3] rantau, 
                                        //          [4] mu, 
                                        //          [5] InvQ, 
                                        //          [6] InvSigma
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
                double* predictor,                         
                int* K,           // IN [1] number of classes
                int* N,           // IN [1] number of observations
                int* n,           // IN [1] number of subjects
                int* nY,          // IN [3]        counts of Nums, Ords and Bins variables
                int* nfix,        // IN [sum(nY)]  number of FIXED  regressors for each response
                int* cumnfix,     // IN [totnY] cummulative number of fixed  regressors
                int* totnfix,     // IN [1] total number of fixed  regressors
                int* n_j,         // IN [n] number of observations dedicated to j-th subject
                int* i_j,         // IN [n * max_n_j] indeces dedicated to j-th subject 
                int* nk           // IN [K] number of observations in classes 1, ..., K
                )    
{ 
  /** Declarations **/
  int j, k, l;                  // looping indeces
  int y;                        // variable index
  int coord;                    // coordinate for matrices as arrays
  double newa, newb;            // updated set of parameters for gamma distribution
  double x;                     // auxiliary value for squares
  
  // Is tau class-specific?
  if(spec[1]){
    // also spec["beta"]=T
    
    for(k = 0; k < *K; k++){
      for(y = 0; y < nY[0]; y++){
        // prior contribution
        newa = *((*param).a);
        newb = *((*param).b);
        // model contributions
        newa += (nk[k] + nfix[y])/2.0;
        // beta prior contribution
        for(j = 0; j < nfix[y]; j++){
          x = (*last).beta[j + cumnfix[y] + k*dims[7]] - (*param).fixmu[j + cumnfix[y]];
          newb += 0.5 * x * x/(*param).fixD[j + cumnfix[y]];
        }
        // model contribution - only those observations in k-th class
        for(j = 0; j < *n; j++){
          if((*last).U[j] == k){
            for(l = 0; l < n_j[j]; l++){
              coord = i_j[j + *n * l];
              x = Y[coord + *N * y] - predictor[coord + *N * y];
              newb += 0.5 * x * x;
            }
          }
        }
        
        // generating new rantau
        (*last).tau[y + k * dims[8]] = rgamma(newa, 1/newb);
      }
    }

  }else{ 
    // there is only one tau for each numeric response to be generated
    
    for(y = 0; y < nY[0]; y++){
      // prior contribution
      newa = *((*param).a) + *N/2.0;
      newb = *((*param).b);
      // beta prior contributions
      if(spec[2]){
        // beta is class-specific
        newa += *K * nfix[y]/2.0;
        
        for(k = 0; k < *K; k++){
          for(j = 0; j < nfix[y]; j++){
            x = (*last).beta[j + cumnfix[y] + k*dims[7]] - (*param).fixmu[j + cumnfix[y]];
            newb += 0.5 * x * x/(*param).fixD[j + cumnfix[y]];
          }
        }
        
      }else{
        // beta is not class-specific
        newa += nfix[y]/2.0;
        
        for(j = 0; j < nfix[y]; j++){
          x = (*last).beta[j + cumnfix[y]] - (*param).fixmu[j + cumnfix[y]];
          newb += 0.5 * x * x/(*param).fixD[j + cumnfix[y]];
        }
      }
      
      // model contribution - all observations
      for(j = 0; j < *N; j++){
          x = Y[j + *N * y] - predictor[j + *N * y];
          newb += 0.5 * x * x;
      }
      
      // generating new tau
      (*last).tau[y] = rgamma(newa, 1/newb);
      //printf("y=%d, newa = %f, newb = %f\n", y, newa, newb);
    }
    
  } // end of if spec["tau"]
  
  
} // end of void 
