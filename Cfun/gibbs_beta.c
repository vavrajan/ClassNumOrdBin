/*
 * Generating betas from the conditioned distribution 
 */


#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "matrmult.h"
#include "cholesky.h"

#include "structures.h"

void gibbs_beta(struct str_state* last,  // OUT+IN last known values of generated parameters
                struct str_param* param, // IN hyperparameters
                int* Id,              // [N] IDs
                double* Y,            // IN [*N * (1 + totnY)]
                double* X,            // [N*(1 + #regr)]    ID + regressors
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
                int* N,           // IN [1] number of observations
                int* n,           // IN [1] number of subjects
                int* nY,          // IN [3]        counts of Nums, Ords and Bins variables
                int* totnY,       // IN [1] total number of responses
                int* nfix,        // IN [sum(nY)]  number of FIXED  regressors for each response
                int* nran,        // IN [sum(nY)]  number of RANDOM regressors for each response
                int* cumnfix,     // IN [totnY] cummulative number of fixed  regressors
                int* cumnran,     // IN [totnY] cummulative number of random regressors
                int* totnfix,     // IN [1] total number of fixed  regressors
                int* totnran,     // IN [1] total number of random regressors
                int* FormulaF,    // IN [sum(nY)]  numbers of columns of X that should be used for FIXED  effects of modelled responses
                int* FormulaR,    // IN [sum(nY)]  numbers of columns of X that should be used for RANDOM effects of modelled responses
                int* n_j,         // IN [n] number of observations dedicated to j-th subject
                int* i_j,         // IN [n * max_n_j] indeces dedicated to j-th subject 
                int* nUk          // IN [K] number of subjects     in classes 1, ..., K
                )    
{ 
  /** Declarations **/
  int j, k, l, m, p, q;             // looping indeces
  int id;                           // id of the subject
  int y;                            // variable index
  int coord;                        // coordinate for matrices as arrays
  int row, col, col1, col2;         // indeces to orient in the symmetric matrix
  double chol[*totnfix * (*totnfix + 1)/2];        // cholesky decomposition of matrix M
  //double betahat[*totnfix];         // mean value of the normal distribution
  double almost_betahat[*totnfix];  // chol^{-T} * bhat --> to be backsolved to get muhat
  double rightb[*totnfix];          // right-hand side of the equation to get muhat
  double rvec[*totnfix];            // vector of generated N(0,1) values
  double scaled_rvec[*totnfix];     // vector of generated N(0,Sigma) values
  double XtX[*totnfix * (*totnfix+1)/2]; // XtX + diag(1/param.fixD)
  double modY[*N];                  // modified response (Y/latent - b*Z)
  double* pY;                       // pointer to Y
  double tau;                       // tau value to multiply XtX
  
  for(y = 0; y < *totnY; y++){
    // pointer to the y-th variable
    if(y < nY[0]){
      // numeric variable
      pY = Y + *N * y;
    }else{
      // ordinal or binary variable
      pY = (*last).latent + *N * (y - nY[0]);
    }
    
    // modified response = response - rand effects* regressors
    for(j = 0; j < *N; j++){
      // cycle through all observations
      id = Id[j]; // converting the type double into int
      modY[j] = pY[j]; 
      
      for(m = 0; m < nran[y]; m++){
        modY[j] -= (*last).b[id + *n * (cumnran[y] + m)] * X[j + *N * (FormulaR[cumnran[y] + m])];
      }
    }
    
    
        
    // Is beta class-specific?
    if(spec[2]){
      for(k = 0; k < *K; k++){
        // Is tau class specific? 
        if(y < nY[0]){ // numeric variable
          if(spec[1]){
            tau = (*last).tau[y + k * dims[8]];
          }else{
            tau = (*last).tau[y];
          }
        }else{ // ordinal or binary variable
          tau = 1;
        }
        // XtX (for class k) + D calculation
        for(p = 0; p < nfix[y]; p++){
          for(q = p; q < nfix[y]; q++){
            coord = p + q*(q+1)/2;
            // start with fixD and then add XtX
            XtX[coord] = (p == q) / (*param).fixD[p+cumnfix[y]]; // 0 off diagonal, 1/fixD on the diagonal
            
            for(j = 0; j < *n; j++){
              // cycle through subjects... but only those that are in the class k
              if((*last).U[j] == k){
                for(l = 0; l < n_j[j]; l++){
                  // cycle through observations dedicated to subject j
                  row = i_j[j + *n * l];
                  col1 = FormulaF[cumnfix[y] + p];
                  col2 = FormulaF[cumnfix[y] + q];
                  XtX[coord] += X[row + *N * col1] * X[row + *N * col2];
                }
              }
            }
            XtX[coord] *= tau;
          }
        } // end of XtX + 1/fixD calculation
        
        // rightb calculation: Rcode: t(pomXfix)%*%tildeY[selection[[k]]]+param$fixmu[[y]]/param$fixD[[y]]
        for(p = 0; p < nfix[y]; p++){
          rightb[p] = (*param).fixmu[p+cumnfix[y]] / (*param).fixD[p+cumnfix[y]];
          // now add t(X) * modY
          for(j = 0; j < *n; j++){
            // cycle through subjects... but only those that are in the class k
            if((*last).U[j] == k){
              for(l = 0; l < n_j[j]; l++){
                // cycle through observations dedicated to subject j
                row = i_j[j + *n * l];
                col = FormulaF[cumnfix[y] + p];
                rightb[p] += X[row + *N * col] * modY[row];
              }
            }
          }
          rightb[p] *= tau;
        }

        // mean value calculation
        cholesky_solve(XtX, nfix + y, rightb, chol, almost_betahat);
        //backsolve2(chol, almost_betahat, nfix + y, betahat);
        
        // now we have mean value (betahat) and also the inverse of variance matrix (XtX)
        for(l = 0; l < nfix[y]; l++){
          //rvec[l] = rnorm(0.0, 1.0); // generate from N(0,1) distribution
          rvec[l] = almost_betahat[l] + rnorm(0.0, 1.0);
        }
        // scale it 
        backsolve2(chol, rvec, nfix + y, scaled_rvec);
        
        // add mean value and store to last
        for(l = 0; l < nfix[y]; l++){
          //(*last).beta[l + k * dims[7] + cumnfix[y]] = betahat[l] + scaled_rvec[l];
          (*last).beta[l + k * dims[7] + cumnfix[y]] = scaled_rvec[l];
          //printf("%f, ", scaled_rvec[l]);
        }
        
      } // end of for k
    }else{
      // beta is not class specific --> tau cannot be as well
      if(y < nY[0]){ // numeric variable
        tau = (*last).tau[y];
      }else{ // ordinal or binary variable
        tau = 1;
      }
      // last$beta[[y]] <- backsolve(cholXtXplusD[[y]], rnorm(nfix[y])) + as.numeric(solve(XtXplusD[[y]], t(Xfix[[y]])%*%tildeY+param$fixmu[[y]]/param$fixD[[y]]))
      // if(whatsave["beta"]){chain$beta[[y]][i,] <- last$beta[[y]]}
      
      // XtX (for class k) + D calculation
      for(p = 0; p < nfix[y]; p++){
        for(q = p; q < nfix[y]; q++){
          coord = p + q*(q+1)/2;
          // start with fixD and then add XtX
          XtX[coord] = (p == q) / (*param).fixD[p+cumnfix[y]]; // 0 off diagonal, 1/fixD on the diagonal
          
          for(j = 0; j < *n; j++){
            // cycle through subjects... but only those that are in the class k
            for(l = 0; l < n_j[j]; l++){
              // cycle through observations dedicated to subject j
              row = i_j[j + *n * l];
              col1 = FormulaF[cumnfix[y] + p];
              col2 = FormulaF[cumnfix[y] + q];
              XtX[coord] += X[row + *N * col1] * X[row + *N * col2];
            }
          }
          XtX[coord] *= tau;
        }
      } // end of XtX + 1/fixD calculation
      
      // rightb calculation
      for(p = 0; p < nfix[y]; p++){
        rightb[p] = (*param).fixmu[p+cumnfix[y]] / (*param).fixD[p+cumnfix[y]];
        // now add t(X) * modY
        for(j = 0; j < *n; j++){
          // cycle through subjects... but only those that are in the class k
          for(l = 0; l < n_j[j]; l++){
            // cycle through observations dedicated to subject j
            row = i_j[j + *n * l];
            col = FormulaF[cumnfix[y] + p];
            rightb[p] += X[row + *N * col] * modY[row];
          }
        }
        rightb[p] *= tau;
      }
      
      
      // mean value calculation
      cholesky_solve(XtX, nfix + y, rightb, chol, almost_betahat);
      //backsolve2(chol, almost_betahat, nfix + y, betahat);
      
      // now we have mean value (betahat) and also the inverse of variance matrix (XtX)
      for(l = 0; l < nfix[y]; l++){
        //rvec[l] = rnorm(0.0, 1.0); // generate from N(0,1) distribution
        rvec[l] = almost_betahat[l] + rnorm(0.0, 1.0);
      }
      // scale it 
      backsolve2(chol, rvec, nfix + y, scaled_rvec);
      
      // add mean value and store to last
      for(l = 0; l < nfix[y]; l++){
        //(*last).beta[l + cumnfix[y]] = betahat[l] + scaled_rvec[l];
        (*last).beta[l + cumnfix[y]] = scaled_rvec[l];
      }
    } // end of if(spec["beta"])
    
  } // end of for y
  
  
} // end of void 
