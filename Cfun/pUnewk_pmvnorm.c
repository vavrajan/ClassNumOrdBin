//
//  PURPOSE:   Calculation of classification probabilitites for 
//              subject not originally in the study
//
//  AUTHOR:    Jan Vavra
//             vavraj[AT]karlin.mff.cuni.cz
//
//  LOG:       20191031  created
//
// ================================================================
//
// To obtain a dll file in Windows, run
//
//    library("callr")
//    rcmd(cmd = "SHLIB", cmdargs = "pUnewk_pmvnorm.c")
//
//    To be precise, you need to compile all in one:
// rcmd(cmd = "SHLIB", cmdargs = c("pUnewk_calculation.c",
//                                  "calculate_single_predictor.c",
//                                  "cholesky.c",
//                                  "matrmult.c",
//                                  "mypmvnorm.c))
//
// To use declared functions in R use
//
//    dyn.load("./pUnewk_pmvnorm.dll")
//
  
  
  /*** These are headers available within the R source tree      ***/
  /*** that provide mathematical and also many statistical       ***/
  /*** functions (densities, cdf's, random number generators)    ***/
/*** available in R itself.                                    ***/

#include <R.h>
#include <Rmath.h>
#include <Rinternals.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>
#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <mvtnorm.h>
#include <mvtnormAPI.h>

#include "calculate_single_predictor.h"
#include "matrmult.h"
#include "cholesky.h"
#include "mypmvnorm.h"

#define LARGE_VALUE 1000000000000
#define LOW_VALUE  -1000000000000
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))


/*** ================================================================================ ***/
/*** THE KEY PART OF THE CODE ***/
/*** Generating of b_new and latent variables ***/
/*** And then the classification probabilities ***/
/*** ================================================================================ ***/

void pUnewk_pmvnorm( 
    /** IN parameters **/
    int* Id,              // [Nnew] IDs
    double* Y,            // [Nnew*sum(nY)]  responses
    double* X,            // [N*#regr]    regressors
    int* spec,            // [7]                class-specific parameters
    // order:   [0] gamma, 
    //          [1] tau, 
    //          [2] beta, 
    //          [3] rantau, 
    //          [4] mu, 
    //          [5] InvQ, 
    //          [6] InvSigma
    /** Parameters describing dimensions **/
    int* chain,           // [1]        number of generated chain (just for printing)
    int* K,               // [1]        number of classes
    int* BM,              // [1]        total number of generated states
    int* FormulaF,        // [sum(nY)]  numbers of columns of X that should be used for FIXED  effects of modelled responses
    int* FormulaR,        // [sum(nY)]  numbers of columns of X that should be used for RANDOM effects of modelled responses
    int* Nnew,            // [1]        total number of observations
    int* dims,            // [7]       the length of subarray that corresponds to one state (disected by various parameters)
    // order  [0]   w,
    //        [1]   pUnewk,
    //        [2]   gamma
    //        [3]   beta,
    //        [4]   tau,
    //        [5]   mu,
    //        [6]   InvSigma
    int* dimswithK,       // [7]       the length of subarray that corresponds to one state 
    //            (disected by various parameters, also multiplication by K incorporated when such parameters is class-specific)
    // order  [0]   w,
    //        [1]   pUnewk,
    //        [2]   gamma
    //        [3]   beta,
    //        [4]   tau,
    //        [5]   mu,
    //        [6]   InvSigma
    
    int* ncat,            // [nY[1 .. Ords]]  the counts of categories of ordinal variables
    int* nY,              // [3]        counts of Nums, Ords and Bins variables
    int* nnew,            // [1]        total number of subjects (different ids in the dataset)
    int* nfix,            // [sum(nY)]  number of FIXED  regressors for each response
    int* nran,            // [sum(nY)]  number of RANDOM regressors for each response
    //double* predictor,
    // those are used to construct totnfix, totnran and cummulative versions
    /** Arrays to store generated states **/
    /** Some of them might be NULL (see whatsave and calc) **/
    double* w,            // [(B+M) * K]               Class probabilities
    double* pUnewk,       // [(B+M) * nnew*K]          Probability of i-th subject belong to each class
    double* pUerror,      // [(B+M) * nnew*K]          Estimated absolute error, with 99% confidence level
    double* pUinform,     // [(B+M) * nnew*K]          Occured error when calculating probability using pmvnorm
                              //  if INFORM = 0, normal completion with ERROR < EPS;
                              //  if INFORM = 1, completion with ERROR > EPS and MAXPTS function values used; 
                              //                 increase MAXPTS to decrease ERROR;
                              //  if INFORM = 2, N > 1000 or N < 1.
                              //  if INFORM = 3, correlation matrix not positive semi-definite.
    double* gamma,        // [(B+M) * sum(ncat-2)]     Division points 
    double* gamma1,       // [nY[1]]  fixed thresholds for     Ordinal variables
    double* gammaBin,     // [1]      fixed threshold  for all Binary  variables
    double* beta,         // [(B+M) * totnfix]         Beta parameters for fixed effects
    double* tau,          // [(B+M) * nY[0]]           Precision parameter of numerical variables
    double* mu,           // [(B+M) * totnran]         Mean values of random effects
    double* InvSigma)     // [(B+M) * totnran*(totnran+1)/2] Precision matrix of random effects
    
{
  /*** Declarations ***/
  /*** ------------ ***/
  int i, j, l, k, m;    // looping indeces
  int y;                // index variable for outcomes
  
  /*** Calculation of useful dimensions, etc. ***/
  /*** -------------------------------------- ***/
  // Total number of responses
  int totnY = 0;            // total number of responses
  for(i = 0; i < 3; i++)
    totnY += nY[i];
  
  // Total and cummulative numbers of fixed and random regressors
  int totnfix = 0;      // total number of fixed regressors
  int totnran = 0;      // total number of random regressors
  int cumnfix [totnY + 1];  // cummulative number of fixed and random regressors
  int cumnran [totnY + 1];  // cummulative number of fixed and random regressors
  
  cumnfix[0] = cumnran[0] = 0;
  
  for(y = 0; y < totnY; y++){
    totnfix += nfix[y];
    totnran += nran[y];
    cumnfix[y+1] = totnfix;
    cumnran[y+1] = totnran;
  }
  
  // Bounds and gamma parameter
  int maxnbounds = 0; // maximal number of needed bounds for truncating latent ordinal variables
  for(y = 0; y < nY[1]; y++){
    if(maxnbounds < ncat[y]-1){
      maxnbounds = ncat[y]-1;
    }
  }
  
  // Subjects 
  int n_j [*nnew];      // number of observations dedicated to j-th subject
  int max_n_j = 0;      // maximal value of n_j numbers
  int id;               // ID of the subject
  
  for(i = 0; i < *nnew; i++){
    n_j[i] = 0;
  }
  
  for(i = 0; i < *Nnew; i++){
    id = Id[i]; // (conversion of type double into int)
    n_j[id]++; 
    if(n_j[id] > max_n_j){
      max_n_j = n_j[id];
    }
  }
  
  
  int i_j [*nnew * max_n_j];      // matrix of indeces to j-th subject
  int nn_jj [*nnew];              // counts of each subject
  
  for(i = 0; i < *nnew; i++){
    nn_jj[i] = 0;
  }
  for(i = 0; i < *Nnew; i++){
    id = Id[i];
    i_j[id + *nnew * nn_jj[id]] = i;
    ++nn_jj[id];
  }
  
  /*** Main calculation cycle ***/
  /***   i = current state    ***/
  /*** ---------------------- ***/
  
  /* Declarations */
  double p [*K];              // probabilities
  double sump;                // sum of probabilities
  double logp;                // logarithm of some probability
  
  int percentage = 0;         // How many percent has been already generated?
  int newperc = 0;
  
  double* pw;           // pointer to w
  double* pgamma;       // pointer to gamma
  double* ppbeta;        // pointer to beta
  double* ptau;         // pointer to tau
  double* pmu;          // pointer to mu
  //double* pInvSigma;    // pointer to InvSigma
  
  double chol [dimswithK[6]];         // for Cholesky decomposition of InvSigma
  double* pchol;                      // pointer for chol decomposition
  double cholinv [dimswithK[6]];      // place to store inversions of chol
  double* pcholinv;                   // pointer for cholinv
  double ZC [max_n_j*totnY*dims[5]];  // matrix Z*cholinv (by columns)
  double Mean [max_n_j*totnY];  // vector of means
  double V [max_n_j*totnY * (max_n_j*totnY+1)/2]; // upper triangle of variance matrix
  
  double ynum [max_n_j*nY[0]];        // numeric part of measured outcomes
  double chold1 [max_n_j*nY[0] * (max_n_j*nY[0]+1)/2];            // for Cholesky decomposition of V[1:d1,1:d1]
  double Cy [max_n_j*nY[0]];          // for solution to C^T * x = (Ynum - Mean[1:d1])
  double detVd1;                      // for determinant of the first block of V
  double InvVy [max_n_j*totnY];       // for C^{-1}*C^{-T}*(Y-mean)
  double CV12 [max_n_j*max_n_j*nY[0]*(nY[1]+nY[2])]; // for C{-T}V12
  
  double condmean [max_n_j*(nY[1]+nY[2])]; 
                  // conditional mean of 2nd block given the first block
  double condvar [max_n_j*(nY[1]+nY[2]) * (max_n_j*(nY[1]+nY[2])+1)/2]; 
                  // conditional variance of 2nd block given the first block
  double corr [max_n_j*(nY[1]+nY[2]) * (max_n_j*(nY[1]+nY[2])-1)/2];
                  // correlation matrix of conditional variance of 2nd block
                  // !! for the purpose of pmvnorm needs to be lower tringle by columns!!
  double chold2 [max_n_j*(nY[1]+nY[2]) * (max_n_j*(nY[1]+nY[2])-1)/2];  // for Cholesky decomposition of V[d1:d1+d2,d1:d1+d2]
                  
  
  int cat;                // auxiliary value to store current category in {0,1,2,ncat[y]-1}
  int dimgammay;          // dimension of current part of gamma parameter
  int cumdimgammay;       // cummulative sum of dimensions of gamma parameter
  int infin [max_n_j * (nY[1] + nY[2])]; // type of infinity bounds for integration
  double lower [max_n_j * (nY[1] + nY[2])]; // lower bounds for integration
  double upper [max_n_j * (nY[1] + nY[2])]; // upper bounds for integration
  
  
  int d1, d2;                   // dimensions of block dedicated tu Nums and OrdsBins 
  int pos, bdim, xdim, coord;   // auxiliary dimension identifiers
  
  // parameters for pmvnorm function
  int nu = 0;             // degrees of freedom for student distribution (if 0 then MV Normal)
  int maxpts = 25000;     // default in mvtnorm: 25000
  double releps = 0;      // default in mvtnorm: 0
  double abseps = 0.001;  // default in mvtnorm: 0.001
  double error;           // output value of error 
  double value;           // output value of the integral
  int rnd = 1;            // Get/PutRNGstate (TRUE or FALSE)
  int inform = 0;         // information about pmvnorm process (see explanation above)
  double alpha = 0.01241933; // leads to quantile 2.5 as in pmvnorm
  
  
  double delta[max_n_j * totnY];
  
  for(i = 0; i < (max_n_j * totnY); ++i){
    delta[i] = 0.0;
  }
  
  
  
  for(i = 0; i < *BM; i++){
    // Printing current iteration number
    //printf("Generating state   %d   out of   %d  ...  %d percent \n", i+1, *BM, 100*(i+1)/(*BM));
    newperc = 100*(i+1)/(*BM);
    
    if(percentage < newperc){
      percentage = newperc;
      //printf("Chain: %d, %d\n", *chain, percentage);
      printf("\r");
      printf("Chain: %d, %c%d %c%c", *chain, '(', percentage, '%', ')');
      fflush(stdout);
    }
    
    // pointer actualisation
    // w
    pw = w + i*dims[0];
    
    // cholesky decomposition of InvSigmas 
    // (ahead to not to compute it several times)
    pchol = chol;
    pcholinv = cholinv;
    
    if(spec[6]){
      for(k = 0; k < *K; k++){
        justcholesky(InvSigma + i*dimswithK[6] + k*dims[6], 
                     pchol, dims + 5);
        // inverse calculation
        invofchol(pchol, pcholinv, dims + 5);
        pchol += dims[6];
        pcholinv += dims[6];
      }
    }else{
      justcholesky(InvSigma + i*dimswithK[6], pchol, dims + 5);
      // inverse calculation
      invofchol(pchol, pcholinv, dims + 5);
    }
    
    //if((i == 0) && (j == 4) && (*chain == 1)){
    //  printf("\nChol 1");
    //  for(l = 0; l < dims[5]; l++){
    //    printf("\n");
    //    for(m = 0; m < dims[5]; m++){
    //      if(m<l){
    //        printf("%f, ", 0.0);
    //      }else{
    //        printf("%f, ", chol[l + m*(m+1)/2]);
    //      }
    //    }
    //  }
    //}
    
    
    //if((i == 0) && (j == 4) && (*chain == 1)){
    //  printf("\nInv of chol 1");
    //    for(l = 0; l < dims[5]; l++){
    //    printf("\n");
    //    for(m = 0; m < dims[5]; m++){
    //      if(m<l){
    //        printf("%f, ", 0.0);
    //      }else{
    //        printf("%f, ", cholinv[l + m*(m+1)/2]);
    //      }
    //    }
    //  }
    //}
    
    //printf("\nChol 2");
    //if(i == 0){
    //  for(l = 0; l < dims[5]; l++){
    //    printf("\n");
    //    for(m = 0; m < dims[5]; m++){
    //      if(m<l){
    //        printf("%f, ", 0.0);
    //      }else{
    //        printf("%f, ", cholinv[dims[6] + l + m*(m+1)/2]);
    //      }
    //    }
    //  }
    //}
    
    //printf("\nInv of chol 2");
    //if(i == 0){
    //  for(l = 0; l < dims[5]; l++){
    //    printf("\n");
    //    for(m = 0; m < dims[5]; m++){
    //      if(m<l){
    //        printf("%f, ", 0.0);
    //      }else{
    //        printf("%f, ", cholinv[dims[6] + l + m*(m+1)/2]);
    //      }
    //    }
    //  }
    //}
    
    // Cycle through all subjects
    for(j = 0; j < *nnew; j++){
      // Cycle through all classes
      
      sump = 0;
      for(k = 0; k < *K; k++){
        // zero step - pointer actualisation
        // gamma
        if(spec[0]){
          pgamma = gamma + i*dimswithK[2] + k*dims[2];
        }else{
          pgamma = gamma + i*dims[2];
        }
        // beta
        if(spec[2]){
          ppbeta = beta + i*dimswithK[3] + k*dims[3];
        }else{
          ppbeta = beta + i*dims[3];
        }
        // tau
        if(spec[1]){
          ptau = tau + i*dimswithK[4] + k*dims[4];
        }else{
          ptau = tau + i*dims[4];
        }
        // mu
        if(spec[4]){
          pmu = mu + i*dimswithK[5] + k*dims[5];
        }else{
          pmu = mu + i*dims[5];
        }
        // InvSigma
        if(spec[6]){
          //pInvSigma = InvSigma + i*dimswithK[6] + k*dims[6];
          pchol = chol + k*dims[6];
          pcholinv = cholinv + k*dims[6];
        }else{
          //pInvSigma = InvSigma + i*dims[6];
          pchol = chol;
          pcholinv = cholinv;
        }
        
        //// Step 1 - start with w
        p[k] = pw[k];
        
        //// Preparations for Step 2
        // Block dimensions
        d1 = n_j[j] * nY[0];          // dim of block Y_Nums
        d2 = n_j[j] * (nY[1]+nY[2]);  // dim of block Y_OrdsBins
        
        // computing Mean of combined continuous variables
        for(y = 0; y < totnY; y++){
          for(l = 0; l < n_j[j]; l++){
            pos = y*n_j[j] + l;
            Mean[pos] = 0;
            // add fixed part
            for(m = 0; m < nfix[y]; m++){
              bdim = m + cumnfix[y];
              xdim = i_j[j + *nnew * l]  + *Nnew * FormulaF[bdim];
              // add m-th regressor for y multiplied by corresponding beta
              Mean[pos] += X[xdim] * ppbeta[bdim];
            }
            // add random part
            for(m = 0; m < nran[y]; m++){
              bdim = m + cumnran[y];
              xdim = i_j[j + *nnew * l]  + *Nnew * FormulaR[bdim];
              // add m-th regressor for y multiplied by corresponding mu
              Mean[pos] += X[xdim] * pmu[bdim];
            }
          }
        }
        //if((i == 0) && (j == 4) && (*chain == 1)){
        //  printf("\n Mean = ");
        //  for(y = 0; y < totnY; y++){
        //    for(l = 0; l < n_j[j]; l++){
        //      printf("%f, ", Mean[y*n_j[j] + l]);
        //    }
        //  }
        //}
        
        // computing V of combined continuous variables
        // inverse of cholesky decomposition of InvSigma already calculated
        // calculate Z*C first and then (Z*C) * (Z*C)^T
        for(y = 0; y < totnY; y++){
          for(l = 0; l < n_j[j]; l++){
            for(m = 0; m < dims[5]; m++){
              coord = l + y*n_j[j] + m*n_j[j]*totnY;
              ZC[coord] = 0;
              for(int ii = cumnran[y]; (ii < (cumnran[y] + nran[y])) && (ii <= m); ii++){
                // ii is restricted by the upper-right triangularity of cholinv
                xdim = i_j[j + *nnew * l] + *Nnew * FormulaR[ii];
                ZC[coord] += X[xdim] * pcholinv[ii + m*(m+1)/2];
              }
            }
          }
        }
        
        //if((i == 0) && (j == 4) && (*chain == 1)){
        //  printf("\n ZC");
        //  for(l = 0; l < (d1+d2); l++){
        //    printf("\n");
        //   for(m = 0; m < dims[5]; m++){
        //      printf("%f, ", ZC[l + m*(d1+d2)]);
        //    }
        //  }
        //}
        
        // and now (Z*C) * (Z*C)^T
        for(int y1 = 0; y1 < totnY; y1++){
          for(int l1 = 0; l1 < n_j[j]; l1++){
            for(int y2 = y1; y2 < totnY; y2++){
              for(int l2 = ((y2>y1)?(0):(l1)); l2 < n_j[j]; l2++){
                // couple (y,l) sets one coordinate
                coord = l1 + y1*n_j[j] + (l2 + y2*n_j[j])*(l2 + y2*n_j[j]+1)/2;
                V[coord] = 0;
                for(m = MAX(cumnran[y1], cumnran[y2]); m < dims[5]; m++){
                  // not necessary to count from m=0 when all before that max would be zeros
                  V[coord] += ZC[l1 + y1*n_j[j] + m*n_j[j]*totnY] * ZC[l2 + y2*n_j[j] + m*n_j[j]*totnY];
                }
              }
            }
          }
        }
        // dont forget on the diagonal
        for(y = 0; y < totnY; y++){
          for(l = 0; l < n_j[j]; l++){
            coord = (l + y*n_j[j])*(3 + l + y*n_j[j])/2;
            if(y < nY[0]){
              V[coord] += 1/ptau[y];  // numeric variable
            }else{
              V[coord] += 1;          // Ord or Bin variable
            }
          }
        }
        
        //if((i==0) && (j==4) && (*chain == 1)){
        //  printf("\n V");
        //  for(l = 0; l < (d1); l++){
        //    printf("\n");
        //    for(m = 0; m < (d1); m++){
        //     if(m < l){
        //        printf("%f, ", 0.0);
        //      }else{
        //        printf("%f, ", V[l + m*(m+1)/2]);
        //      }
        //    }
        //  }
        //}
        
        // Now Mean and V are the parameters of normal distribution
        // Y_Nums and Y_latent |params    (random effects b integrated out)
          
        //// Step 2 - multiply by p(Y_Nums|params)
        // in R: dmvnorm(yy[1:d1], mean = Mean[1:d1], sigma = V[1:d1,1:d1])
        for(y = 0; y < nY[0]; y++){
          for(l = 0; l < n_j[j]; l++){
            ynum[l + y*n_j[j]] = Y[i_j[j + *nnew * l] + *Nnew * y] - Mean[l + y*n_j[j]];
          }
        }
        // first block of Mean and V is in the front d1 or d1*(d1+1)/2 values
        //cholesky2(V, &d1, ynum, chold1, Cy, &detVd1); // wrong, does C^{-1}*ynum, not C^{-T}*ynum
        
        //if((i == 0) && (j == 4) && (*chain == 1)){
        //  printf("\n ynum   ");
        //  for(l = 0; l < d1; l++){
        //    printf("%f, ", ynum[l]);
        //  }
        //}
        
        justcholesky(V, chold1, &d1);
        
        //if((i==0) && (j==4) && (*chain == 1)){
        //  printf("\n chold1");
        //  for(l = 0; l < (d1); l++){
        //    printf("\n");
        //    for(m = 0; m < (d1); m++){
        //      if(m < l){
        //        printf("%f, ", 0.0);
        //      }else{
        //        printf("%f, ", chold1[l + m*(m+1)/2]);
        //      }
        //    }
        //  }
        //}
        
        forwardsolve2(chold1, ynum, &d1, Cy); // finds solution to C^T * Cy = ynum
        
        //if((i == 0) && (j == 4) && (*chain == 1)){
        //  printf("\n Cy   ");
        //  for(l = 0; l < d1; l++){
        //    printf("%f, ", Cy[l]);
        //  }
        //}
        
        detVd1 = 1;
        for(l = 0; l < d1; l++){
          detVd1 *= chold1[l*(l+3)/2]; // product of elements on the diagonal
        }
        detVd1 *= detVd1;     // square it
        logp = 0;
        for(l = 0; l < d1; l++){
          logp += Cy[l] * Cy[l];
        }
        //if((i==0) && (j==4) && (*chain == 1)){
        //  printf("\nk = %d, ()TV-1() = %f", k, logp);
        //}
        logp += d1*log(2 * M_PI) + log(detVd1); // M_PI is C definition of pi from <math.h>
        // part d1*log(2*pi) could be omitted because of the scaling 
        // p/sum(p) in which it is the same for all p[k]
        
        p[k] *= exp(-0.5*logp);
        
        //if((i==0) && (j==4) && (*chain == 1)){
        //  printf("\nk = %d, p(k) = %f, dmvnorm = %f", k, p[k], exp(-0.5*logp));
        //}
        
        //// Preparations for Step 3
        // now we need to find parameters of conditional distribution of
        //  Y_OrdBin | Y_Num
        /// Conditional mean
        // To the mean we need to add V21 V11^{-1} * (Y-Nums - Mean[1:d1])
        // we already have Cy, what remains to be done is: V21 * chold1^{-1} * Cy
        backsolve2(chold1, Cy, &d1, InvVy); // InvVy = chold1^{-1} * Cy
        
        for(l = 0; l < d2; l++){
          // Initialize as the unconditioned mean value
          condmean[l] = Mean[d1+l];
          // now add V21 * InvVy
          for(m = 0; m < d1; m++){
            // we need V21[l,m] = V12[m,l] element
            condmean[l] += V[m + (l+d1)*(l+d1+1)/2] * InvVy[m];
          }
        }
        
        //if((i == 0) && (j == 4) && (*chain == 1)){
        //  printf("\n condmean = ");
        //  for(l = 0; l < d2; l++){
        //    printf("%f, ", condmean[l]);
        //  }
        //}
        
        /// Conditional variance
        // V22 - V21 * V11^{-1} * V12
        // first start with chold1^{-T}*V12 - column by column
        for(l = 0; l < d2; l++){
          forwardsolve2(chold1, V+(l+d1)*(l+d1+1)/2, &d1, CV12+l*d1);
        }
        
        //if(i < 2 && ((j == 1)|| (j == 6))){
        //  printf("\n CV12");
        //  for(l = 0; l < d1; l++){
        //    printf("\n");
        //    for(m = 0; m < d2; m++){
        //      printf("%f, ", CV12[l + m*d1]);
        //    }
        //  }
        //}
        
        // now compute V22 - CV12^T * CV12
        for(l = 0; l < d2; l++){
          for(m = l; m < d2; m++){
            // initialize by V22 value
            coord = l + m*(m+1)/2;
            condvar[coord] = V[l+d1 + (m+d1)*(m+d1+1)/2];
            // now subtract CV12[column l]^T CV12[column m]
            for(int index = 0; index < d1; index++){
              condvar[coord] -= CV12[index + l*d1] * CV12[index + m*d1];
            }
          }
        }
        
        //if((i == 0) && (j == 4) && (*chain == 1)){
        //  printf("\n condvar");
        //  for(l = 0; l < (d2); l++){
        //   printf("\n");
        //    for(m = 0; m < (d2); m++){
        //      if(m < l){
        //        printf("%f, ", 0.0);
        //      }else{
        //        printf("%f, ", condvar[l + m*(m+1)/2]);
        //      }
        //    }
        //  }
        //}
        
        /// lower and upper bounds calculation
        // infin will declare what type of integration should be performed   
        // INFIN(I) < 0, Ith limits are (-infinity, infinity);
        // INFIN(I) = 0, Ith limits are (-infinity, UPPER(I)];
        // INFIN(I) = 1, Ith limits are [LOWER(I), infinity);
        // INFIN(I) = 2, Ith limits are [LOWER(I), UPPER(I)].
        
        /// Setting bounds + standardization to zero mean and unit variance
        // Ordinal variables first
        
        //if((i == 0) && (j == 2) && (*chain == 1)){
        //  printf("\n j = %d, k = %d, Value = %f, Error = %f, inform = %d",
        //         j, k, value, error, inform);
        //}
        
        cumdimgammay = 0;
        for(y = 0; y < nY[1]; y++){
          dimgammay = ncat[y]-2;
          for(l = 0; l < n_j[j]; l++){
            cat = Y[i_j[j + *nnew * l] + *Nnew * (y+nY[0])]; // in 0,1,2,...,ncat[y]-1
            coord = l + y*n_j[j];
            switch(cat){
              case 0:
                infin[coord] = 0;
                lower[coord] = LOW_VALUE;
                upper[coord] = (gamma1[y] - condmean[coord]);///sqrt(condvar[coord*(coord+3)/2]);
                break;
                
              case 1:
                infin[coord] = 2;
                lower[coord] = (gamma1[y] - condmean[coord]);///sqrt(condvar[coord*(coord+3)/2]);
                upper[coord] = (pgamma[cat-1 + cumdimgammay] - condmean[coord])/sqrt(condvar[coord*(coord+3)/2]);
                break;
                
              default:
                // case with the largest possible category ncat[y]-1 must be modified
                infin[coord] = (((ncat[y]-1) == cat)?(1):(2)); // if last category, then 1, otherwise 2
                lower[coord] = (pgamma[cat-2 + cumdimgammay] - condmean[coord]);///sqrt(condvar[coord*(coord+3)/2]);
                upper[coord] = (((ncat[y]-1) == cat)?(LARGE_VALUE):((pgamma[cat-1 + cumdimgammay] - condmean[coord])));///sqrt(condvar[coord*(coord+3)/2])));
              
            }
            
            //if((i == 0) && (j == 2) && (*chain == 1)){
            //  printf("\n    y = %d, ncat[y] = %d, gamma1[y] = %f, pgamma[0] = %f, pgamma[last] = %f, l = %d, cat = %d, coord = %d, pos = %d, Infin = %d, lower = %f, upper = %f",
            //         y, ncat[y], gamma1[y], pgamma[cumdimgammay], pgamma[cumdimgammay+ncat[y]-3], l, cat, coord, pos, infin[coord], lower[coord], upper[coord]);
            //}
            
            //if(cat == (ncat[y]-1)){
            //  infin[coord] = 1;
            //  upper[coord] = LARGE_VALUE;
            //}
          }
          cumdimgammay += dimgammay;
        }
        // Binary variables next
        for(y = 0; y < nY[2]; y++){
          for(l = 0; l < n_j[j]; l++){
            cat = Y[i_j[j + *nnew * l] + *Nnew * (y+nY[0]+nY[1])];
            coord = l + (y+nY[1])*n_j[j];
            switch(cat){
              case 0:
                infin[coord] = 0;
                lower[coord] = LOW_VALUE;
                upper[coord] = (*gammaBin - condmean[coord]);///sqrt(condvar[coord*(coord+3)/2]);
                break;
                
              case 1:
                infin[coord] = 1;
                lower[coord] = (*gammaBin - condmean[coord]);///sqrt(condvar[coord*(coord+3)/2]);
                upper[coord] = LARGE_VALUE;
                break;
                
            }
            
            //if((i == 0) && (j == 2) && (*chain == 1)){
            //  printf("\n    y = %d, ncat[y] = %d, l = %d, cat = %d, coord = %d, pos = %d, Infin = %d, lower = %f, upper = %f",
            //         y+nY[1], 2, l, cat, coord, pos, infin[coord], lower[coord], upper[coord]);
            //}
          }
        }
        
        // correlation matrix calculation
        // !!dont forget that pmvnorm requires lower triangle stored by columns!!
        for(l = 0; l < (d2-1); l++){
          for(m = l+1; m < d2; m++){
            corr[l+(m-1)*m/2] = condvar[l + m*(m+1)/2]/sqrt(condvar[l*(l+3)/2]*condvar[m*(m+3)/2]);
          }
        }

        //// Step 3 - multiply by pmvnorm(Ystar_OrdsBin|Y_Nums,params)
        // mvtnorm_C_mvtdst is defined in mvtnorm/inst/include/mvtnormAPI.h 
        mvtnorm_C_mvtdst(&d2, &nu, lower, upper,
                         infin, corr, delta,
                         &maxpts, &abseps, &releps,
                         &error, &value, &inform, &rnd);
        
        //justcholesky(condvar, chold2, &d2);                 
        //mypmvnorm(&d2, lower, upper, infin,
        //          chol, &maxpts, &abseps, &alpha,
        //          &error,&value);
        p[k] *= value;
        
        //if((i == 0) && (j == 2) && (*chain == 1)){
        //  printf("\n value = %f, p[k] = %f", value, p[k]);
        //}
        // store error and inform
        pUerror[ i*dims[1] + j + *nnew * k] = error;
        pUinform[i*dims[1] + j + *nnew * k] = inform;
        
        // add p[k] to the sum of probabilities
        sump += p[k];
        
        //if((i == 0) && (j == 4) && (*chain == 1)){
        //  printf("\n j = %d, k = %d, Value = %f, Error = %f, inform = %d",
        //       j, k, value, error, inform);
        //  for(y = 0; y < (nY[1]+nY[2]); y++){
        //    for(l = 0; l < n_j[j]; l++){
        //      printf("\n    y = %d, l = %d, cat = %f, Infin = %d, lower = %f, upper = %f",
        //             y, l, Y[i_j[j + *nnew * l] + *Nnew * (y+nY[0])], infin[l + y*n_j[j]], lower[l + y*n_j[j]], upper[l + y*n_j[j]]);
        //    }
        //  }
        //}

      } // end of for k
      
      for(k = 0; k < *K; k++){
        pUnewk[i*dims[1] + j + *nnew * k] = p[k] / sump;
      }
      
      
      
    } // end of for j
 
  } // end of for i 
  
  printf("\n");
  fflush(stdout);
  
}
