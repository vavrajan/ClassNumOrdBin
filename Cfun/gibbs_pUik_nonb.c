/*
 * Generating probabilities of i-th subject being in class k (random effects b integrated out)
 * Calculated only if calc[0] == TRUE
 */


#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>

#include "matrmult.h"
#include "cholesky.h"

#include "structures.h"

#define LOGPINF -10000000
#define SHIFTLOGP 10

void gibbs_pUik_nonb( struct str_state* last,  // OUT+IN last known values of generated parameters
                      struct str_param* param, // IN hyperparameters
                      double* Y,            // IN [*N * totnY]
                      double* X,            // [N*(#regr)]    ID + regressors
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
                      int* n,           // IN [1] number of subjects
                      int* K,           // IN [1] number of classes
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
                      int* isnonzero,   // IN+OUT [n*K] 0/1 values whether probability is non-zero
                      double* logp_base,// IN [n*K] base for logp, might have been calculated in past
                      int* nUk,         // IN [K] number of subjects     in classes 1, ..., K
                      int* nk,          // IN [K] number of observations in classes 1, ..., K
                      double* pUik_nonb,// IN [(B+M)*dims[14]] the place where to store probabilities
                      double* ZtZ,      // IN [n * dimZtZ] all ZtZ matrices by subjects and responses
                      int* dimZtZ       // IN [1] dimension of one block of all matrices corresponding to one subject
                      )    
{ 
  /** Declarations **/
  int j, k, l, m;         // looping indeces
  int isnonzeroj [*K];    // j-th subjects isnonzero values
  int y;                  // y-th response
  int I;                  // category of l-th observation
  double bound;           // lower bound for categorization of latent ordinal variable
  int sumcat = 0;         // sum of categorical values
  double logp[*K];        // probabilities on logarithmic scale  
  double* pbeta;          // pointer to the begginning of vector beta
  double* ptau;           // pointer to the begginning of vector tau
  double* pInvSigma;      // pointer to the begginning of matrix InvSigma
  double* pmu;            // pointer to the begginning of vector mu
  double pomtau;          // auxiliary value of tau
  double x;               // double value for saving the computation time
  double maxlogp;         // maximal value of logp
  double sumexplogp;      // sum exponentials of logarithms of probabilities
  double explogp[*K];     // exponential of logarithms of probabilities
  double logp_base_forallj [*K];    // logp_base part that is common to all j
  double modY [*max_n_j * *totnY];  // modificated (latent) response
  double matrsum [*totnran * (*totnran + 1)/2]; // sum of two inverse variance matrices
  double vecZY [*totnran];// vector Z times modY ... and then added some other vector
  int row, col;           // row or column indeces - for working with symmetric matrices
  //int index;              // auxiliary index for block modification of matrsum
  int dimnrany;           // (cummulative) dimension of blocks of ZtZ by response 
  double bhat [*totnran]; // vector of mean value of full conditioned distribution 
  double chol [*totnran * (*totnran+1)/2]; // matrix for cholesky decomposition 
  double det;             // determinant of a matrix being Cholesky decomposed
  
  /* Common part of logp_base to all j */
 
  //printf("dimZtZ = %d\n", *dimZtZ);
  //printf("i = %d\n", *i);
  //printf("totnran = %d\n", *totnran);
  //printf("totnY = %d\n", *totnY);
 
  for(k = 0; k < *K; k++){
    // w contribution
      logp_base_forallj[k] = log((*last).w[k]);
    
    // InvSigma contribution (if class-specific)  
      if(spec[6]){
        // InvSigma is class-specific
        logp_base_forallj[k] += 0.5*log((*last).detInvSigma[k]);
      }
    }
  
  /* j-specific part */
  
  for(j = 0; j < *n; j++){
    
    maxlogp = LOGPINF;
    
    for(k = 0; k < *K; k++){
      isnonzeroj[k] = 1;
      // if gamma depends on k, then there might be some zero issues
      if(spec[0]){
        // now iterate only through ordinal variables
        for(y = 0; y < nY[1]; y++){
          for(l = 0; l < n_j[j]; l++){
            I = 1;
            bound = (*param).gamma1[y];
            
            while(bound < (*last).latent[i_j[j + *n * l] + *N * y] 
                    && 
                      I < ncat[y]){
              bound = (*last).gamma[(I-1) + sumcat + k * dims[3]];
              I++;
            }
            // I now has the value of category which belongs to this latent value of y
            // It must be the same as the measured value of y
            // Otherwise the probability will be zero 
            
            isnonzeroj[k] *= Y[i_j[j + *n * l] + *N * y] == I;
            if (isnonzeroj[k] == 0){
              break; // when zero, it does not make sence to continue
            }
          }
          if (isnonzeroj[k] == 0){
            break; // when zero, it does not make sence to continue
          }
        sumcat += ncat[y] - 2;
        }
        if (isnonzeroj[k] == 0){
          break; // when zero, it does not make sence to continue
        }
      } // end of else of if spec["gamma"]
    
    // Store isnonzero output for later use in pUik  
    isnonzero[j + *n * k] = isnonzeroj[k];
      
    // logp_base calculation
    logp_base[j + *n * k] = logp_base_forallj[k]; 
    
    // tau contribution (if tau class-specific - depends on subject j)
    if(spec[1]){
      // tau is class-specific
      for(l = 0; l < dims[8]; l++){
        logp_base[j + *n * k] += 0.5*n_j[j]*log((*last).tau[l + k*dims[8]]);
      }
    }
    
    if(!isnonzeroj[k]){
      logp_base[j + *n * k] = LOGPINF;
    }
    
   
    // logp calculation
    if(isnonzeroj[k]){
      logp[k] = logp_base[j + *n * k];
      
      /* Pointer actualisation */
      // pointer to beta
      pbeta = (*last).beta;
      // Is beta class-specific?
      if(spec[2]){
        pbeta += k*dims[7];
      }
      
      // pointer to tau;
      ptau= (*last).tau;
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
          pomtau = ptau[y];
        }else{
          pomtau = 1.0;
        }
            
        /* matrsum update - needs to be done by blocks */
        for(col = 0; col < nran[y]; col++){
          for(row = 0; row <= col; row++){
            /*if(j==1){
              printf("col = %d, row = %d\n", col, row);
              printf("matrsum dim = %d\n", row + cumnran[y] + (col + cumnran[y])*(col + cumnran[y]+1)/2);
              printf("ZtZ dim = %d\n", *dimZtZ * j + dimnrany + row + col*(col+1)/2);
            }*/
            matrsum[row + cumnran[y] + (col + cumnran[y])*(col + cumnran[y]+1)/2] += pomtau * ZtZ[*dimZtZ * j + dimnrany + row + col*(col+1)/2];
          }
          //index += col + 1 + cumnran[y];
        }

        /* vecZY update - is done by y-blocks */
        for(l = 0; l < nran[y]; l++){
          // to the l-th member add scalar product of column of X and modY (for j-th subject)
          for(m = 0; m < n_j[j]; m++){
            vecZY[l + cumnran[y]] += X[i_j[j + *n * m] + *N * FormulaR[cumnran[y] + l]] * modY[m + y * n_j[j]];
          }
          // multiply by tau
          vecZY[l + cumnran[y]] *= pomtau;
        }
        /* increase auxiliary dimension counters */
        //index += nran[y];
        dimnrany += nran[y] * (nran[y]+1)/2;
        
      } // end of for y (all responses)
       
      /* bhat calculation */
      // bhat solves equation matrsum * bhat = vecZY
      cholesky2(matrsum, totnran, vecZY, chol, bhat, &det);
      
      /* Now we can add specific things to logp */
      // t(bhat) * matrsum * bhat
      for(l = 0; l < *totnran; l++){
        logp[k] += 0.5 * bhat[l]*vecZY[l];
      }
      //aBa(bhat, matrsum, &x, totnran);
      //logp[k] += 0.5 * x;
      
      // tau contribution within modY * modY
      if(spec[1]){
        // tau is class-specific
        for(y = 0; y < *totnY; y++){
          // current tau value
          if(y < nY[0]){
            pomtau = *(ptau + y);
          }else{
            pomtau = 1;
          }
          
          for(l = 0; l < i_j[j]; l++){
            logp[k] -= 0.5 * pomtau * modY[l + y * n_j[j]] * modY[l + y * n_j[j]];
          }
        }
      }
      
      if(spec[4] || spec[6]){
        aBa(pmu, pInvSigma, &x, totnran);
        logp[k] -= 0.5 * x;
        }
      
      // determinant from integration
      logp[k] -= 0.5 * log(det);
      // det = determinant of a symmetric matrix --> NEEDS TO BE IMPLEMENTED or found

    }else{
      logp[k] = LOGPINF;
    } // end of if isnonzeroj[k]
   
    if(maxlogp < logp[k]){
      maxlogp = logp[k];
    }
      
    } // end of for k --> we now have maxlogp
   
    sumexplogp = 0;
    for(k = 0; k < *K; k++){
      logp[k] += SHIFTLOGP - maxlogp;
      explogp[k] = exp(logp[k]);
      sumexplogp += explogp[k];
    }
    
    for(k = 0; k < *K; k++){
      explogp[k] /= sumexplogp;
      pUik_nonb[*i * dims[14] + j + *n * k] = explogp[k];
    }
   
  } // end of for j  
    

  
} // end of void 
