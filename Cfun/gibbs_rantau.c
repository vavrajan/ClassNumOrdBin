/*
 * Generating rantau from the conditioned distribution 
 */


#include <R.h>
#include <Rmath.h>

#include "structures.h"

void gibbs_rantau(struct str_state* last,  // OUT+IN last known values of generated parameters
                  struct str_param* param, // IN hyperparameters
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
                  int* K,           // IN [1] number of classes
                  int* totnY,       // IN [1] total number of responses
                  int* nran,        // IN [sum(nY)]  number of FIXED  regressors for each response
                  int* cumnran,     // IN [totnY] cummulative number of fixed  regressors
                  int* totnran      // IN [1] total number of fixed  regressors
                  )    
{ 
  /** Declarations **/
  int k, l;                  // looping indeces
  double newa, newb;            // updated set of parameters for gamma distribution
  double x;                     // auxiliary value for squares
  
  // Is rantau class-specific?
  if(spec[3]){
    // also spec["mu"]=T
    
    for(k = 0; k < *K; k++){
      // prior contribution
      newa = *((*param).aran) + 0.5;
      
      // prior of mu contribution to newb + generating
      for(l = 0; l < *totnran; l++){
        newb = *((*param).bran);
        x = (*last).mu[l + k*dims[9]] - (*param).ranmu0[l];
        newb += 0.5 * x * x;
        
        // generating new rantau
        (*last).tau[l + k*dims[9]] = rgamma(newa, 1/newb);
      } // end of for l
    } // end of for k

  }else{ 
    // rantau is not class-specific, however, mu still can be
    
    for(l = 0; l < *totnran; l++){
      newa = *((*param).aran);
      newb = *((*param).bran);
      // is mu class specific?
      if(spec[4]){
        // mu is class-specific
        // prior mu contribution to new a
        newa += *K/2.0;
        // prior mu contribution to new b
        for(k = 0; k < *K; k++){
          x = (*last).mu[l + k*dims[9]] - (*param).ranmu0[l];
          newb += 0.5 * x * x;
        }
      }else{
        // mu is not class-specific
        // prior mu contribution to new a
        newa += 0.5;
        // prior mu contribution to new b
        x = (*last).mu[l] - (*param).ranmu0[l];
        newb += 0.5 * x * x;
      }
      // now we have newa and new b for l-th parameter of y-th variable

      // generating new rantau
      (*last).rantau[l] = rgamma(newa, 1/newb);
    } // end of for l
    
  } // end of if spec["rantau"]
  
  
} // end of void 
