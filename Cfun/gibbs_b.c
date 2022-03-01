/*
 * Generating random effects from the full conditioned distribution 
 */


#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "matrmult.h"
#include "cholesky.h"

#include "structures.h"


void gibbs_b( struct str_state* last,  // OUT+IN last known values of generated parameters
              struct str_param* param, // IN hyperparameters
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
              int* FormulaF,    // IN [sum(nY)]  numbers of columns of X that should be used for FIXED  effects of modelled responses
              int* FormulaR,    // IN [sum(nY)]  numbers of columns of X that should be used for RANDOM effects of modelled responses
              int* N,           // IN [1] number of observations 
              int* n,           // IN  [1] number of subjects
              int* K,           // IN  [1] number of classes
              int* ncat,        // IN [nY[1 .. Ords]]  the counts of categories of ordinal variables
              int* nY,          // IN [3]        counts of Nums, Ords and Bins variables
              int* totnY,       // IN [1] total number of responses
              int* nfix,        // IN [sum(nY)]  number of FIXED  regressors for each response
              int* nran,        // IN [sum(nY)]  number of RANDOM regressors for each response
              int* cumnfix,     // IN [totnY] cummulative number of fixed  regressors
              int* cumnran,     // IN [totnY] cummulative number of random regressors
              int* totnfix,     // IN [1] total number of fixed  regressors
              int* totnran,     // IN [1] total number of random regressors
              int* i,           // IN [1] iteration number
              int* n_j,         // IN [n] number of observations dedicated to j-th subject
              int* max_n_j,     // IN [1] maximal number of all n_j values
              int* i_j,         // IN [n * max_n_j] indeces dedicated to j-th subject 
              double* ZtZ,      // IN [n * dimZtZ] all ZtZ matrices by subjects and responses
              int* dimZtZ       // IN [1] dimension of one block of all matrices corresponding to one subject
              )    
{ 
  /** Declarations **/
  int j, k, l, m;         // looping indeces
  int y;                  // variable index
  double* pbeta;          // pointer to the begginning of vector beta
  double* ptau;           // pointer to the begginning of vector tau
  double* pInvSigma;      // pointer to the begginning of matrix InvSigma
  double* pmu;            // pointer to the begginning of vector mu
  double pomtau;          // auxiliary value of tau
  double v [*totnran];    // vector of rnorm(0,1) values 
  double out [*totnran];  // vector of rnorm(0,Sigma) values 
  double modY [*max_n_j * *totnY];  // modificated (latent) response
  double matrsum [*totnran * (*totnran + 1)/2]; // sum of two inverse variance matrices
  double vecZY [*totnran];// vector Z times modY ... and then added some other vector
  int row, col;           // row or column indeces - for working with symmetric matrices
  //int index;              // auxiliary index for block modification of matrsum
  int dimnrany;           // (cummulative) dimension of blocks of ZtZ by response 
  //double bhat [*totnran]; // vector of mean value of full conditioned distribution 
  double chol [*totnran * (*totnran+1)/2]; // matrix for cholesky decomposition 
  //double det;             // determinant of matrix being Cholesky decomposed
  
  /* We will calculate for each subject j separately */
  
  for(j = 0; j < *n; j++){
    
    k = (*last).U[j];
    
    /* Pointer actualisation */
    // pointer to beta
    pbeta = (*last).beta;
    // Is beta class-specific?
    if(spec[2]){
      pbeta += k*dims[7];
    }
    
    // pointer to tau;
    ptau = (*last).tau;
    if(spec[1]){
      ptau += k * dims[8];
    }
        
    // pointer to InvSigma
    pInvSigma = (*last).InvSigma;
    if(spec[6]){
      pInvSigma += k*dims[11];
    }
    
    // pointer to mu
    pmu = (*last).mu;
    if(spec[4]){
      pmu += k*dims[9];
    }
    
    /* matrsum inicialization */
    for(l = 0; l < dims[11]; l++){
      matrsum[l] = pInvSigma[l];
    }
    
    /* vecZY inicialization */
    for(l = 0; l < *totnran; l++){
      vecZY[l] = 0;
      for(m = 0; m < *totnran; m++){
        // row = min(m,l), col = max(m,l)
        if(m < l){
          row = m;
          col = l;
        }else{
          row = l;
          col = m;
        }
        vecZY[l] += pInvSigma[row + col*(col+1)/2] * pmu[m];
      }
    }
      
    //index = 0;
    dimnrany = 0;
          
    for(y = 0; y < *totnY; y++){
      
      for(l = 0; l < n_j[j]; l++){
        // modified value of Y = its numeric version - X^T beta
        if(y < nY[0]){
          // it is a numerical variable
          modY[l + y * n_j[j]] = Y[i_j[j + *n * l] + *N * y];
        }else{
          // it is an ordinal or binary variable --> use latent variable instead
          modY[l + y * n_j[j]] = (*last).latent[i_j[j + *n * l] + *N * (y-nY[0])];
        }

        // subtract scalar product of X and beta
        for(m = 0; m < nfix[y]; m++){
          modY[l + y * n_j[j]] -= X[i_j[j + *n * l] + *N * FormulaF[cumnfix[y] + m]] * pbeta[cumnfix[y] + m];
        }
      } // end of for l (all observations with j-th subject)
      
      /* tau value for such y */
      if(y < nY[0]){
        pomtau = *(ptau + y);
      }else{
        pomtau = 1;
      }
          
      /* matrsum update - needs to be done by blocks */
      for(col = 0; col < nran[y]; col++){
        for(row = 0; row <= col; row++){
          //matrsum[index + row] += pomtau * ZtZ[*dimZtZ * j + dimnrany + row + col*(col+1)/2];
          matrsum[cumnran[y] + row + (cumnran[y]+col)*(cumnran[y]+col+1)/2] += pomtau * ZtZ[*dimZtZ * j + dimnrany + row + col*(col+1)/2];
        }
        //index += col + 1 + cumnran[y];
      }
      
      /* vecZY update - is done by y-blocks */
      for(l = 0; l < nran[y]; l++){
        // to the l-th member add scalar product of column of X and modY (for j-th subject)
        for(m = 0; m < n_j[j]; m++){
          vecZY[l + cumnran[y]] += pomtau * X[i_j[j + *n * m] + *N * FormulaR[cumnran[y] + l]] * modY[m + y * n_j[j]];
        }
        // multiply by tau
        //vecZY[l + cumnran[y]] *= pomtau;
      }

      /* increase auxiliary dimension counters */
      //index += nran[y];
      dimnrany += nran[y] * (nran[y]+1)/2;
      
    } // end of for y (all responses)
      
    /* bhat calculation */
    // bhat solves equation matrsum * bhat = vecZY
    cholesky_solve(matrsum, totnran, vecZY, chol, v);
    //backsolve2(chol, v, totnran, bhat);
    
    /* Generating random effects */

    for(l = 0; l < *totnran; l++){
      //v[l] = rnorm(0.0, 1.0);
      v[l] += rnorm(0.0, 1.0);
    }
    
    // scaling
    backsolve2(chol, v, totnran, out);
    
    // shifting
    for(l = 0; l < *totnran; l++){
      //(*last).b[j + *n * l] = out[l] + bhat[l];
      (*last).b[j + *n * l] = out[l];
    }
    
  } // end of for j  
    

  
} // end of void 
