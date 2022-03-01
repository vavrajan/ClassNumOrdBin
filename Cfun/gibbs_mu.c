/*
 * Generating mu from the conditioned distribution 
 */


#include <R.h>
#include <Rmath.h>

#include "matrmult.h"
#include "cholesky.h"

#include "structures.h"

void gibbs_mu(struct str_state* last,  // OUT+IN last known values of generated parameters
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
              int* K,                 // IN [1] number of classes
              int* n,                 // IN [1] number of subjects
              int* totnran,           // IN [1] total number of random regressors
              int* nUk                // IN [K] number of subjects     in classes 1, ..., K
              )    
{ 
  /** Declarations **/
  int j, k, l;                  // looping indeces
  int coord;                    // coordinate for matrices as arrays
  int row, col;                 // indeces to orient in the symmetric matrix
  double Smatrix[dims[11]];     // Scale matrix of Wishart distribution
  double chol[dims[11]];        // cholesky decomposition of matrix M
  double* pInvSigma;            // pointer to InvSigma  
  double* prantau;              // pointer to rantau 
  //double muhat[dims[9]];        // mean value of the normal distribution
  double almost_muhat[dims[9]]; // chol^{-T} * bhat --> to be backsolved to get muhat
  double rightb[dims[9]];       // right-hand side of the equation to get muhat
  double bhat[dims[9]];         // mean value of the random effects
  double rvec[dims[9]];         // vector of generated N(0,1) values
  double scaled_rvec[dims[9]];  // vector of generated N(0,Sigma) values
  
  // Is mu class-specific?
  if(spec[4]){
    // there is *K different values of mu to be generated
    
    for(k = 0; k < *K; k++){
      /* Pointer settings */
      // InvSigma
      pInvSigma = (*last).InvSigma;
      if(spec[6]){
        pInvSigma += k*dims[11];
      }
      // rantau
      prantau = (*last).rantau;
      if(spec[3]){
        prantau += k*dims[10];
      }
      
      // combination of variance matrices
      for(l = 0; l < *totnran; l++){
        // diagonal
        Smatrix[l*(l+3)/2] = prantau[l] + nUk[k] * pInvSigma[l*(l+3)/2];
        // off-diagonal
        for(j = l+1; j < *totnran; j++){
          coord = l + j*(j+1)/2;
          Smatrix[coord] = nUk[k] * pInvSigma[coord];
        }
      }
      
      // bhat calculation - average of last random effects b
      for(l = 0; l < dims[9]; l++){
        bhat[l] = 0;
        for(j = 0; j < *n; j++){
          // cycle through all subjects, but choose only those in k-th class
          if((*last).U[j] == k){
            bhat[l] += (*last).b[j + *n * l];
          }
        }
        bhat[l] /= nUk[k]; // nUk should contain the total number of subjects being of the value k
      }
      
      // rightb calculation (right-hand side)
      for(j = 0; j < dims[9]; j++){
        rightb[j] = prantau[j] * (*param).ranmu0[j];
        // now the InvSigma * bhat part
        for(l = 0; l < dims[9]; l++){
          if(l < j){
            row = l; // min
            col = j; // max
          }else{
            row = j; // min
            col = l; // max
          }
          rightb[j] += nUk[k] * pInvSigma[row + col*(col + 1)/2] * bhat[l];
        }
      }
      
      // muhat is the solution to Smatrix * x = rightb
      cholesky_solve(Smatrix, dims + 9, rightb, chol, almost_muhat);
      //backsolve2(chol, almost_muhat, dims + 9, muhat);
      
      // now we have mean value (muhat) and also the inverse of variance matrix (Smatrix)
      for(l = 0; l < dims[9]; l++){
        //rvec[l] = rnorm(0.0, 1.0); // generate from N(0,1) distribution
        rvec[l] = almost_muhat[l] + rnorm(0.0, 1.0);
      }
      // scale it 
      backsolve2(chol, rvec, dims + 9, scaled_rvec);
      
      // add mean value and store to last
      for(l = 0; l < dims[9]; l++){
        //(*last).mu[l + k * dims[9]] = muhat[l] + scaled_rvec[l];
        (*last).mu[l + k * dims[9]] = scaled_rvec[l];
      }

    } // end of for k
    
  }else{ // then rantau cannot be class-specific
    // there is only one mu to be generated
    // rantau
    prantau = (*last).rantau;
    
    // Is InvSigma class-specific?
    if(spec[6]){
      // Smatrix is the weighted sum of InvSigma matrices and diag(rantau)
      for(l = 0; l < *totnran; l++){
        for(j = l; j < *totnran; j++){
          coord = l + j*(j+1)/2;
          Smatrix[coord] = (l == j) * prantau[l];
          for(k = 0; k < *totnran; k++){
            Smatrix[coord] += nUk[k] * (*last).InvSigma[coord + k*dims[11]];
          }
        }
      }
      
      // bhat calculation - average of last random effects b
      for(k = 0; k < *K; k++){
        for(l = 0; l < dims[9]; l++){
          coord = l + k*dims[9];
          bhat[coord] = 0;
          for(j = 0; j < *n; j++){
            // cycle through all subjects, but choose only those in k-th class
            if((*last).U[j] == k){
              bhat[coord] += (*last).b[j + *n * l];
            }
          }
          bhat[coord] /= nUk[k]; // nUk should contain the total number of subjects being of the value k
        }
      }
      
      // rightb calculation (right-hand side)
      for(j = 0; j < dims[9]; j++){
        rightb[j] = prantau[j] * (*param).ranmu0[j];
        // now the InvSigma * bhat part
        for(k = 0; k < dims[9]; k++){
          for(l = 0; l < dims[9]; l++){
            if(l < j){
              row = l; // min
              col = j; // max
            }else{
              row = j; // min
              col = l; // max
            }
            rightb[j] += nUk[k] * (*last).InvSigma[row + col*(col + 1)/2 + k*dims[11]] * bhat[l + k * dims[9]];
          }
        }
      }
      
    }else{
      // InvSigma is not class-specific
      // Smatrix is the weighted sum of InvSigma matrix and diag(rantau)
      for(l = 0; l < *totnran; l++){
        Smatrix[l*(l+3)/2] = prantau[l] + *n * (*last).InvSigma[l*(l+3)/2];
        for(j = l+1; j < *totnran; j++){
          coord = l + j*(j+1)/2;
          Smatrix[coord] = *n * (*last).InvSigma[coord];
        }
      }
      
      // bhat calculation - average of last random effects b
      for(l = 0; l < dims[9]; l++){
        bhat[l] = 0;
        for(j = 0; j < *n; j++){
          // cycle through all subjects -- we do not care about the class
          bhat[l] += (*last).b[j + *n * l];
        }
        bhat[l] /= *n; // divide by the total number of subjects
      }
      
      // rightb calculation (right-hand side)
      for(j = 0; j < dims[9]; j++){
        rightb[j] = prantau[j] * (*param).ranmu0[j];
        // now the InvSigma * bhat part
        for(l = 0; l < dims[9]; l++){
          if(l < j){
            row = l; // min
            col = j; // max
          }else{
            row = j; // min
            col = l; // max
          }
          rightb[j] += *n * (*last).InvSigma[row + col*(col + 1)/2] * bhat[l];
        }
      }
    } // end of if spec["Invsigma"]
    // the rest goes the same way regardless of spec["InvSigma"]
    
    // muhat is the solution to Smatrix * x = rightb
    cholesky_solve(Smatrix, dims + 9, rightb, chol, almost_muhat);
    //backsolve2(chol, almost_muhat, dims + 9, muhat);
    
    // now we have mean value (muhat) and also the inverse of variance matrix (Smatrix)
    for(l = 0; l < dims[9]; l++){
      rvec[l] = almost_muhat[l] + rnorm(0.0, 1.0); 
      //rvec[l] = rnorm(0.0, 1.0); // generate from N(0,1) distribution
    }
    // scale it 
    backsolve2(chol, rvec, dims + 9, scaled_rvec);
    
    // add mean value and store to last
    for(l = 0; l < dims[9]; l++){
      //(*last).mu[l] = muhat[l] + scaled_rvec[l];
      (*last).mu[l] = scaled_rvec[l];
    }
    
  } // end of if spec["mu"]
  
  
} // end of void 
