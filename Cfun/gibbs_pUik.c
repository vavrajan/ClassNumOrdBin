/*
 * Generating probabilities of i-th subject being in class k
 * and also the next categories U
 */


#include <R.h>
#include <Rmath.h>
#include <math.h>

#include "matrmult.h"
#include "structures.h"

#define LOGPINF -10000000
#define SHIFTLOGP 10

void gibbs_pUik(struct str_state* last,  // OUT+IN last known values of generated parameters
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
                int* nk           // IN [K] number of observations in classes 1, ..., K
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
  double pompred [*max_n_j]; // X*beta values of specific y and subject j
  double x;               // double value for saving the computation time
  double v [*totnran];    // vector of double values 
  double maxlogp;         // maximal value of logp
  double sumexplogp;      // sum exponentials of logarithms of probabilities
  double explogp[*K];     //exponential of logarithms of probabilities
  int newU;               // newly generated U
  double upom;            // uniformly generated value from interval 0, 1
  double sump;            // cummulative sum of probabilities
  int was_zero_before;   // 0/1 value
  
  // to be updated
  for(k = 0; k < *K; k++){
    nUk[k] = 0;
    nk[k]  = 0;
  }
  
  // cycle through all subjects j
  for(j = 0; j < *n; j++){
    
    if(*i > 0 && calc[0]){
      for(k = 0; k < *K; k++){
        isnonzeroj[k] = isnonzero[j + *n * k];
        }
    }else{
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
              
              isnonzeroj[k] *= (Y[i_j[j + *n * l] + *N * y] == I);
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
      } // end of for k
    } // end of else(i>1 & calc["pUik_nonb"])
    
    for(k = 0; k < *K; k++){
      isnonzero[j + *n * k] = isnonzeroj[k];
    }
    
    // now we know which parts of probs should be positive and which should be zero
    // if isnonzeroj[k]==F --> probability of U_j==k is zero
    if(*i > 0 && calc[0]){ 
      // logp_base was calculated during calculation of in previous iteration pUik_nonb
      for(k = 0; k < *K; k++){
        logp[k] = logp_base[j + *n * k];
      }
    }else{
      for(k = 0; k < *K; k++){
        if(isnonzeroj[k]){
          logp[k] = log((*last).w[k]);
          if(spec[1]){
            // tau is class-specific
            for(l = 0; l < dims[8]; l++){
              logp[k] += 0.5*n_j[j]*log((*last).tau[l + k*dims[8]]);
            }
          }
          if(spec[6]){
            // InvSigma is class-specific
            logp[k] += 0.5*log((*last).detInvSigma[k]);
          }
        }else{
          logp[k] = LOGPINF;
        }
        logp_base[j + *n * k] = logp[k];
      } // end of for k
    } // end of else(i > 0 & calc["pUik_nonb"])
    

    
    // we need to add what needs to be added to the base that is common to pUik and pUik_nonb
    maxlogp = LOGPINF;
    for(k = 0; k < *K; k++){
      if(isnonzeroj[k]){
        // pointer to beta
        pbeta = (*last).beta;
        // is beta class-specific?
        if(spec[2]){
          pbeta += k*dims[7];
        }
        
        /* Model contribution */
        // all numeric, ordinal and binary variables at once
        if(spec[1] || spec[2]){
          for(y = 0; y < *totnY; y++){
            
            // pointer to tau actualisation
            ptau = (*last).tau + y;
            if(spec[1]){
              ptau += k*dims[8];
            }
            // tau is 1 for non-numeric variables
            if(y >= nY[0]){
              pomtau = 1.0;
            }else{
              pomtau = *ptau;
            }
            
            for(l = 0; l < n_j[j]; l++){
              pompred[l] = 0;
              // fixed effects
              for(m = 0; m < nfix[y]; m++){
                pompred[l] += X[i_j[j + *n * l] + *N * FormulaF[cumnfix[y] + m]] * pbeta[cumnfix[y] + m];
              }
              // random effects
              for(m = 0; m < nran[y]; m++){
                pompred[l] += X[i_j[j + *n * l] + *N * FormulaR[cumnran[y] + m]] * (*last).b[j + *n * (cumnran[y] + m)];
              }
              // subtracting the weighted squares
              if(y < nY[0]){
                x = pompred[l] - Y[i_j[j + *n * l] + *N * y];
              }else{
                x = pompred[l] - (*last).latent[i_j[j + *n * l] + *N *(y-nY[0])];
              }
              logp[k] -= 0.5 * pomtau * x * x;
            }
          } // end of for y 
        } // end of if beta or tau are class-specific
        
        /* random effect contribution */
        
        if(spec[4] || spec[6]){
          pInvSigma = (*last).InvSigma;
          if(spec[6]){ 
            pInvSigma += k*dims[11];
          }
          pmu = (*last).mu;
          if(spec[4]){
            pmu += k*dims[9];
          }
          
          for(l = 0; l < *totnran; l++){
            v[l] = pmu[l] - (*last).b[j + *n * l];
          }
          aBa(v, pInvSigma, &x, totnran);
          logp[k] -= 0.5 * x;
        }
        
      }else{ // k-th probability will be zero
        logp[k] = LOGPINF;
      }
      
      // finding the maximum
      if(maxlogp < logp[k]){
        maxlogp = logp[k];
      }
    } // end of for k --> now we have maxlogp, so that we can subtract it
    sumexplogp = 0;
    for(k = 0; k < *K; k++){
      logp[k] += SHIFTLOGP - maxlogp;
      explogp[k] = exp(logp[k]);
      sumexplogp += explogp[k];
    }
    
    for(k = 0; k < *K; k++){
      explogp[k] /= sumexplogp;
      (*last).pUik[j + *n * k] = explogp[k];
    }
    
    // new state U
    newU = 0;
    sump = explogp[0];
    was_zero_before = !isnonzeroj[0];
    upom = runif(0.0, 1.0);
    while(sump < upom){
      newU++;
      sump += explogp[newU];
      was_zero_before = was_zero_before || !isnonzeroj[newU];
    } // newU in {0, 1, ..., K-1}
    
    if(isnonzeroj[newU]){
      (*last).U[j] = newU;
    }else{
      // it still might happen that we want to assign a class, that should not be assigned!!!
      // we will find the class that is closest to the current one
      // it depends wheter we have already met 
      while(!isnonzeroj[newU]){
        if(was_zero_before){ 
          if(newU ==  0) break;
          newU--;
        }else{
          if(newU == *K - 1) break;
          newU++;
        }
      }
      (*last).U[j] = newU;
    }
    
  // update of nUk and nk
  nk[newU] += n_j[j];
  nUk[newU]++;
    
  } // end of for j
  
  
  /*
  // update of parameters
  int id;   // ID number of subject
  
  for(k = 0; k < *K; k++){
    nUk[k] = 0;
    nk[k]  = 0;
  }
  
  for(j = 0; j < *N; j++){
    id = Id[j];        // observation id
    k  = (*last).U[id]; // observations class
    nk[k]++;
  }
  
  for(j = 0; j < *n; j++){
    k = (*last).U[j];  // subjects class
    nUk[k]++;         
  }
   */
  
}
