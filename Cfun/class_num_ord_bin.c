//
//  PURPOSE:   Gibbs algorithm for classification in joint modelling
//             of numeric, ordinal and binary variables
//
//  AUTHOR:    Jan Vavra
//             vavraj[AT]karlin.mff.cuni.cz
//
//  LOG:       20190225  created
//
// ================================================================
//
// To obtain a dll file in Windows, run
//
//    library("callr")
//    rcmd(cmd = "SHLIB", cmdargs = "class_num_ord_bin.c")
//
//    To be precise, you need to compile all in one:
// rcmd(cmd = "SHLIB", cmdargs = c("class_num_ord_bin.c",
//                                  "calculate_predictor.c",
//                                  "cholesky.c",
//                                  "matrmult.c",
//                                  "myrdirichlet.c",
//                                  "myrtruncnorm.c",
//                                  "myrwishart.c",
//                                  "gibbs_b.c",
//                                  "gibbs_beta.c",
//                                  "gibbs_InvQ.c",
//                                  "gibbs_InvSigma.c",
//                                  "gibbs_latent.c",
//                                  "gibbs_mu.c",
//                                  "gibbs_pUik.c",
//                                  "gibbs_pUik_nonb.c",
//                                  "gibbs_rantau.c",
//                                  "gibbs_tau.c"))
//
// To use declared functions in R use
//
//    dyn.load("./class_num_ord_bin.dll")
//


/*** These are headers available within the R source tree      ***/
/*** that provide mathematical and also many statistical       ***/
/*** functions (densities, cdf's, random number generators)    ***/
/*** available in R itself.                                    ***/

#include <R.h>
#include <Rmath.h>
#include <stdio.h>

#include "calculate_predictor.h"
#include "matrmult.h"
#include "cholesky.h"
#include "myrdirichlet.h"

#include "structures.h"

#include "gibbs_b.h"
#include "gibbs_beta.h"
#include "gibbs_InvQ.h"
#include "gibbs_InvSigma.h"
#include "gibbs_latent.h"
#include "gibbs_mu.h"
#include "gibbs_pUik.h"
#include "gibbs_pUik_nonb.h"
#include "gibbs_rantau.h"
#include "gibbs_tau.h"



/*** ================================================================================ ***/
/*** THE KEY PART OF THE CODE ***/
/*** Gibbs algorithm ***/
/*** ================================================================================ ***/

void gibbs_class_num_ord_bin( 
  /** IN parameters **/
  int* Id,              // [N] IDs
  double* Y,            // [N*sum(nY)]  responses
  double* X,            // [N*#regr]    regressors
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
  int* whatsave,        // [14]:              what parameters are saved
                          // order:   [0-2]   w, U, PUik, 
                          //          [3-6]   gamma, min, max, latent, 
                          //          [7-8]   beta, tau, 
                          //          [9-10]  mu, rantau, 
                          //          [11-12] InvSigma, InvQ, 
                          //          [13]    b
  double* vecparam,     // vector of parameters --> str_param
  double* vecinits,     // vector of initial values (except U) --> str_state 
  int* inits_U,         // initial values for U
  /** OUT last state **/
  double* veclast,      // vector of last values (except U) <-- str_state 
  int* last_U,          // last values for U
  /** Parameters describing dimensions **/
  int* chain,           // [1]        number of generated chain (just for printing)
  int* K,               // [1]        number of classes
  int* BM,              // [1]        total number of generated states
  int* FormulaF,        // [sum(nY)]  numbers of columns of X that should be used for FIXED  effects of modelled responses
  int* FormulaR,        // [sum(nY)]  numbers of columns of X that should be used for RANDOM effects of modelled responses
  int* N,               // [1]        total number of observations
  int* dims,            // [23]       the length of subarray that corresponds to one state (disected by various parameters)
                            // order  [0-2]   w, U, pUik,
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
                            // order  [0-2]   w, U, pUik,
                            //        [3-6]   gamma, min, max, latent,
                            //        [7-8]   beta, tau,
                            //        [9-10]  mu, rantau,
                            //        [11-12] InvSigma, InvQ,
                            //        [13]    b,
                            //        [14]    pUik_nonb,
                            //        [15-16] Sigma, Q,
                            //        [17-20] det(Sigma,InvSigma,Q,InvQ),
                            //        [21-22] sdSigma, corSigma,
                            
  int* ncat,            // [nY[1 .. Ords]]  the counts of categories of ordinal variables
  int* nY,              // [3]        counts of Nums, Ords and Bins variables
  int* n,               // [1]        total number of subjects (different ids in the dataset)
  int* nfix,            // [sum(nY)]  number of FIXED  regressors for each response
  int* nran,            // [sum(nY)]  number of RANDOM regressors for each response
  int* dimZtZ,         // [1]        total memory for ZtZ (block diagonal matrices)
  //double* ZtZ,        
  //double* predictor,
  // those are used to construct totnfix, totnran and cummulative versions
  /** Arrays to store generated states **/
  /** Some of them might be NULL (see whatsave and calc) **/
  double* w,            // [(B+M)*K]                    Class probabilities
  int*    U,            // [(B+M)*n]                    To which class does i-th subject belong to
  double* pUik,         // [(B+M)*n * K]                Probability of i-th subject belong to each class
  double* latent,       // [(B+M)*(nY[1]+nY[2])*N]      Latent variables
  double* gamma,        // [(B+M)*sum(ncat-2)(* K)]     Division points 
  double* min,          // [(B+M)*sum(ncat-2)(* K)]     Minimal value after  division point
  double* max,          // [(B+M)*sum(ncat-2)(* K)]     Maximal value before division point 
  double* beta,         // [(B+M)*totnfix(* K)]         Beta parameters for fixed effects
  double* tau,          // [(B+M)*nY[0](* K)]           Precision parameter of numerical variables
  double* mu,           // [(B+M)*totnran(* K)]         Mean values of random effects
  double* rantau,       // [(B+M)*totnran(* K)]         Prior precision parameters of mu
  double* InvSigma,     // [(B+M)*totnran*(totnran+1)/2(* K)] Precision matrix of random effects
  double* InvQ,         // [(B+M)*totnran*(totnran+1)/2(* K)] Prior inverse scale matrix of InvSigma
  double* b,            // [(B+M)*totnran * n]          Random effects
  double* pUik_nonb,    // [(B+M)*n * K]                Probability of i-th subject belong to each class (integrated b)
  double* Sigma,        // [(B+M)*totnran*(totnran+1)/2(* K)] Variance matrix of random effects
  double* Q,            // [(B+M)*totnran*(totnran+1)/2(* K)] Prior scale matrix of InvSigma
  double* detSigma,     // [(B+M)(* K)]                 determinant of Sigma
  double* detInvSigma,  // [(B+M)(* K)]                 determinant of InvSigma
  double* detQ,         // [(B+M)(* K)]                 determinant of Q
  double* detInvQ,      // [(B+M)(* K)]                 determinant of InvQ
  double* sdSigma,      // [(B+M)*totnran(* K)]         standard deviations of random effects
  double* corSigma      // [(B+M)*totnran*(totnran-1)/2(* K)] correlations of random effects
  )
{
  /*** Declarations ***/
  /*** ------------ ***/
  int i, j, l;//k;          // looping indeces
  int row, col;             // looping indeces for matrix dimensions
  
  
  /*** Calculation of useful dimensions, etc. ***/
  /*** -------------------------------------- ***/
  // Total number of responses
  int totnY;            // total number of responses
  totnY = 0;
  for(i = 0; i < 3; i++)
    totnY += nY[i];
  
  // Total and cummulative numbers of fixed and random regressors
  int totnfix = 0;      // total number of fixed regressors
  int totnran = 0;      // total number of random regressors
  int cumnfix [totnY + 1];  // cummulative number of fixed and random regressors
  int cumnran [totnY + 1];  // cummulative number of fixed and random regressors
  
  cumnfix[0] = cumnran[0] = 0;
  
  for(i = 0; i < totnY; i++){
    totnfix += nfix[i];
    totnran += nran[i];
    cumnfix[i+1] = totnfix;
    cumnran[i+1] = totnran;
  }
  
  // Bounds and gamma parameter
  int maxnbounds = 0; // maximal number of needed bounds for truncating latent ordinal variables
  for(int y = 0; y < nY[1]; y++){
    if(maxnbounds < ncat[y]-1){
      maxnbounds = ncat[y]-1;
    }
  }
  
  // Subjects 
  int n_j [*n];       // number of observations dedicated to j-th subject
  int max_n_j = 0;        // maximal value of n_j numbers
  int id;             // ID of the subject
  
  for(i = 0; i < *n; i++){
    n_j[i] = 0;
  }

  for(i = 0; i < *N; i++){
    id = Id[i]; // (conversion of type double into int)
    n_j[id]++; 
    if(n_j[id] > max_n_j){
      max_n_j = n_j[id];
    }
  }
  
  
  int i_j [*n * max_n_j];      // matrix of indeces to j-th subject
  int nn_jj [*n];             // counts of each subject
  
  for(i = 0; i < *n; i++){
    nn_jj[i] = 0;
  }
  for(i = 0; i < *N; i++){
    id = Id[i];
    i_j[id + *n * nn_jj[id]] = i;
    ++nn_jj[id];
  }
  
  /*** Reading parameter values ***/
  /*** ------------------------ ***/
  
  double* p = vecparam; // pointer for elements of vecparam
  struct str_param param;     // structure where hyperparameters will be stored
  
  param.nu_0     = p;     p++;
  param.nu_1     = p;     p++;
  param.delta    = p;     p++;
  param.a        = p;     p++;
  param.b        = p;     p++;
  param.aran     = p;     p++;
  param.bran     = p;     p++;
  param.ranmu0   = p;     p+=totnran;
  param.gamma1   = p;     p+=nY[1];
  param.gammaBin = p;     p++;
  param.fixD     = p;     p+=totnfix;
  param.fixmu    = p;     p+=totnfix;
  param.InvV     = p;     
  
  
  /*** Reading initial values ***/
  /*** ---------------------- ***/
  struct str_state inits;    // pointer to structure where initial values are stored
  double* pp = vecinits;      // pointer for elements of vecinits
  double chol[dims[11]];
  double space_for_det[*K]; // place for inits of detInvSigma
  int k;              // index for class k = 1, ..., K
  
  inits.w        = pp;      pp+=dimswithK[0];
  inits.U        = inits_U;
  inits.pUik     = pp;      pp+=dimswithK[2];
  inits.gamma    = pp;      pp+=dimswithK[3];
  inits.latent   = pp;      pp+=dimswithK[6];
  inits.beta     = pp;      pp+=dimswithK[7];
  inits.tau      = pp;      pp+=dimswithK[8];
  inits.mu       = pp;      pp+=dimswithK[9];
  inits.rantau   = pp;      pp+=dimswithK[10];
  inits.InvSigma = pp;      pp+=dimswithK[11];
  inits.InvQ     = pp;      pp+=dimswithK[12];
  inits.b        = pp;     
  // determinant calculation
  inits.detInvSigma = space_for_det;
  if(spec[6]){
    for(k = 0; k < * K; k++){
      justcholesky(inits.InvSigma + k*dims[11], chol, &totnran);
      detchol(chol, inits.detInvSigma + k, &totnran); 
    }
  }else{
    justcholesky(inits.InvSigma, chol, &totnran);
    detchol(chol, inits.detInvSigma, &totnran);
  }
  
  /*
  printf("%f\n", *(inits.w));
  for(j = 0; j < dimswithK[11]; j++){
    printf("%f\n", inits.InvSigma[j]);
  }
  printf("%f\n", *(inits.detInvSigma));
  printf("%f, %f, %f, %f, %f\n", inits.b[0],
         inits.b[1], inits.b[2], inits.b[3], inits.b[4]);
  printf("%f\n", *(param.gamma1));
  */
  
  


  /** Calculation of initial values of auxiliary parameters **/
  // totnY was calculated above,
  // we need to pass the adress where this value is stored
  // no need in the case of cumnfix and cumnran (arrays - so it is already a pointer)
  
  double predictor[*N * (nY[0] + nY[1] + nY[2])]; // [*N * totnY] current value of predictor
  
  /* predictor */
  calculate_predictor(
    Id, X, spec, predictor,      
    /** Parameters describing dimensions **/
    FormulaF, FormulaR, N, n, dims, &totnY, nfix, nran, cumnfix, cumnran,     
    /** Arrays with necessary parameters from one state **/
    inits.beta, inits.b, inits.U             
  );
  
    

  /* ZtZ values - for each subject and each response separately */
  // MAYBE THE PROBLEM IS RIGHT HERE - CANNOT ALOCATE A MEMORY OF SIZE
  // THAT DEPENDS ON LOCAL VARIABLE!!!
  double ZtZ[*n * *dimZtZ];  // a very long array to store all ZtZ matrices
  int y;            // counter for responses
  int dimnrany = 0;         // cummulative sum of dimensions of symmetric nran[y] matrices
  int ind;                  // auxiliary index --- to not to be computed again and again
  
  for(j = 0; j < *n; j++){
    dimnrany = 0;
    for(y = 0; y < totnY; y++){
      for(col = 0; col < nran[y]; col++){
        for(row = 0; row <= col; row++){
          ind = *dimZtZ * j + dimnrany + row + col * (col+1)/2;
          ZtZ[ind] = 0;
          // scalar product of row-th  and col-th column of Z = X[j's indeces, y's regressors]
          for(l = 0; l < n_j[j]; l++){
            ZtZ[ind] += X[i_j[j + *n * l] + *N * FormulaR[cumnran[y] + row]] * X[i_j[j + *n * l] + *N * FormulaR[cumnran[y] + col]];
          }
        }
      }
      dimnrany += nran[y] * (nran[y]+1)/2;
    }
  }
  
  /* Printing some values for development purposes
  for(j = 0; j < 3; j++){
    printf("Obs of subject %d: %d\n", j, n_j[j]);
    printf("Obs numbers of subject %d: ", j);
    for(l = 0; l < n_j[j]; l++){
      printf("%d, ", i_j[j + *n * l]);
    }
    printf("\n");
  }
   */

  
  /* parameters for classification */
  int nUk [*K];       // number of subjects     in classes 1, ..., K
  int nk  [*K];       // number of observations in classes 1, ..., K
  
  for(k = 0; k < *K; k++){
    nUk[k] = 0;
    nk[k]  = 0;
  }
  
  for(i = 0; i < *N; i++){
    id = Id[i];       // observation id
    k  = inits.U[id]; // observations class
    nk[k]++;
    }
  
  for(i = 0; i < *n; i++){
    k = inits.U[i];  // subjects class
    nUk[k]++;         
  }
  
  
  
  /*** Inicialization of last state ***/
  struct str_state last;              // parameter of last state - will be graduallly updated 
  double ldetInvSigma[dimswithK[18]]; // space for last detInvSigma 
  
  pp = veclast;
  
  /* 
   * last structure will point to the memory given by input veclast and last_U 
   * will be modified gradually during generating states
   * After the generating ends this memory will contain last generated values -->
   * they will be returned back by the function
   */
  
  last.w        = pp;      pp+=dimswithK[0];
  last.U        = last_U;
  last.pUik     = pp;      pp+=dimswithK[2];
  last.gamma    = pp;      pp+=dimswithK[3];
  last.latent   = pp;      pp+=dimswithK[6];
  last.beta     = pp;      pp+=dimswithK[7];
  last.tau      = pp;      pp+=dimswithK[8];
  last.mu       = pp;      pp+=dimswithK[9];
  last.rantau   = pp;      pp+=dimswithK[10];
  last.InvSigma = pp;      pp+=dimswithK[11];
  last.InvQ     = pp;      pp+=dimswithK[12];
  last.b        = pp; 
  last.detInvSigma = ldetInvSigma;
  
  for(i = 0; i < dimswithK[0]; i++){
    last.w[i] = inits.w[i];
  }
  for(i = 0; i < dimswithK[1]; i++){
    last.U[i] = inits.U[i];
  }
  for(i = 0; i < dimswithK[2]; i++){
    last.pUik[i] = inits.pUik[i];
  }
  for(i = 0; i < dimswithK[3]; i++){
    last.gamma[i] = inits.gamma[i];
  }
  for(i = 0; i < dimswithK[6]; i++){
    last.latent[i] = inits.latent[i];
  }
  for(i = 0; i < dimswithK[7]; i++){
    last.beta[i] = inits.beta[i];
  }
  for(i = 0; i < dimswithK[8]; i++){
    last.tau[i] = inits.tau[i];
  }
  for(i = 0; i < dimswithK[9]; i++){
    last.mu[i] = inits.mu[i];
  }
  for(i = 0; i < dimswithK[10]; i++){
    last.rantau[i] = inits.rantau[i];
  }
  for(i = 0; i < dimswithK[11]; i++){
    last.InvSigma[i] = inits.InvSigma[i];
  }
  for(i = 0; i < dimswithK[12]; i++){
    last.InvQ[i] = inits.InvQ[i];
  }
  for(i = 0; i < dimswithK[13]; i++){
    last.b[i] = inits.b[i];
  }
  for(i = 0; i < dimswithK[18]; i++){
    last.detInvSigma[i] = inits.detInvSigma[i];
  }
  
  
  /*
  double lw[dimswithK[0]];
  int    lU[dimswithK[1]];
  double lpUik[dimswithK[2]];
  double lgamma[dimswithK[3]];
  double llatent[dimswithK[6]];
  double lbeta[dimswithK[7]];
  double ltau[dimswithK[8]];
  double lmu[dimswithK[9]];
  double lrantau[dimswithK[10]];
  double lInvSigma[dimswithK[11]];
  double lInvQ[dimswithK[12]];
  double lb[dimswithK[13]]; 
  double ldetInvSigma[dimswithK[18]];
  
  for(i = 0; i < dimswithK[0]; i++){
    lw[i] = inits.w[i];
  }
  for(i = 0; i < dimswithK[1]; i++){
    lU[i] = inits.U[i];
  }
  for(i = 0; i < dimswithK[2]; i++){
    lpUik[i] = inits.pUik[i];
  }
  for(i = 0; i < dimswithK[3]; i++){
    lgamma[i] = inits.gamma[i];
  }
  for(i = 0; i < dimswithK[6]; i++){
    llatent[i] = inits.latent[i];
  }
  for(i = 0; i < dimswithK[7]; i++){
    lbeta[i] = inits.beta[i];
  }
  for(i = 0; i < dimswithK[8]; i++){
    ltau[i] = inits.tau[i];
  }
  for(i = 0; i < dimswithK[9]; i++){
    lmu[i] = inits.mu[i];
  }
  for(i = 0; i < dimswithK[10]; i++){
    lrantau[i] = inits.rantau[i];
  }
  for(i = 0; i < dimswithK[11]; i++){
    lInvSigma[i] = inits.InvSigma[i];
  }
  for(i = 0; i < dimswithK[12]; i++){
    lInvQ[i] = inits.InvQ[i];
  }
  for(i = 0; i < dimswithK[13]; i++){
    lb[i] = inits.b[i];
  }
  for(i = 0; i < dimswithK[18]; i++){
    ldetInvSigma[i] = inits.detInvSigma[i];
  }
  
  last.w = &lw[0];
  last.U = &lU[0];
  last.pUik = &lpUik[0];
  last.gamma = &lgamma[0];
  last.latent = &llatent[0];
  last.beta = &lbeta[0];
  last.tau = &ltau[0];
  last.mu = &lmu[0];
  last.rantau = &lrantau[0];
  last.InvSigma = &lInvSigma[0];
  last.InvQ = &lInvQ[0];
  last.b = &lb[0];
  last.detInvSigma = &ldetInvSigma[0];
  */
  
  /*last = inits; */
  /*
  printf("%f\n", *(last.w));
  for(j = 0; j < dimswithK[11]; j++){
    printf("%f\n", last.InvSigma[j]);
  }
  printf("%f\n", *(last.detInvSigma));
  printf("%f, %f, %f, %f, %f\n", last.b[0],
         last.b[1], last.b[2], last.b[3], last.b[4]);
  */

    
  /*** Gibbs cycle ***/
  /*** ----------- ***/
  
  /* Declarations */
  double alpha [*K];            // vector parameter for Dirichlet distribution
  int isnonzero [*n * *K];      // matrix of 0/1 values for non-zeroability of pUik
  double logp_base [*n * *K];   // logarithm of probabilities - base common to pUik and pUik_nonb
  int percentage = 0;               // How many percent has been already generated?
  int newperc = 0;
  //int sumU = 0;
  
  
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
    
    // probability w
    for(k = 0; k < *K; k++){
      alpha[k] = param.delta[0] + nUk[k];
    }
    // generate new w and store it to last.w
    my_rdirichlet(last.w, K, alpha);
    // saving w
    if(whatsave[0]){
      for(k = 0; k < dims[0]; k++)
        w[i*dims[0] + k] = last.w[k];
    }
    
    // pUik + U
    gibbs_pUik(&last, &param, Y, X, spec, calc, dims, FormulaF, FormulaR, 
               N, n, K, ncat, nY, &totnY, nfix, nran, cumnfix, cumnran, 
               &totnfix, &totnran, &i, n_j, &max_n_j, i_j, isnonzero, logp_base, nUk, nk);
    // i and max_n_j are one-dimensional and NOT pointers
    
    // saving classes U
    if(whatsave[1]){
      for(j = 0; j < dims[1]; j++){
        U[i*dims[1] + j] = last.U[j];
      }
    }
    //sumU=0;
    //for(j = 0; j < dims[1]; j++){
    //  sumU += last.U[j];
    //}
    //printf("sumU = %d: ", sumU);
    //for(j=0; j < 10; j++){
    //  printf("%d, ", last.U[j]);
    //}
    //printf("\n");
    
    // saving probabilities pUik 
    if(whatsave[2]){
      for(j = 0; j < dims[2]; j++){
        pUik[i*dims[2] + j] = last.pUik[j];
      }
    }
    //printf("Puik %d: ", i);
    //for(j=0; j < 10; j++){
    //  printf("%f, ", last.pUik[j]);
    //}
    //printf("\n");
    
    calculate_predictor(
      Id, X, spec, predictor,      
      /** Parameters describing dimensions **/
      FormulaF, FormulaR, N, n, dims, &totnY, nfix, nran, cumnfix, cumnran,     
      /** Arrays with necessary parameters from one state **/
      last.beta, last.b, last.U             
    );
    //printf("Predictor %d: ", i);
    //for(j=0; j < 10; j++){
    //  printf("%f, ", predictor[j]);
    //}
    //printf("\n");
    // Commenting the rest of calculations
    
    // Latent variable modelling parameters (gamma, min, max, latent)
    gibbs_latent(&last, &param, Y, predictor, spec, dims, 
                 N, n, K, ncat, maxnbounds, nY, n_j, i_j, 
                 whatsave + 5, whatsave + 4,
                 max + i*dimswithK[5], min + i*dimswithK[4]);
    if(whatsave[3]){
      for(j = 0; j < dimswithK[3]; j++){
        gamma[i*dimswithK[3] + j] = last.gamma[j];
      }
    }
    //printf("Gamma %d: ", i);
    //for(j=0; j < dimswithK[3]; j++){
    //  printf("%f, ", last.gamma[j]);
    //}
    //printf("\n");
    if(whatsave[6]){
      for(j = 0; j < dimswithK[6]; j++){
        latent[i*dimswithK[6] + j] = last.latent[j];
      }
    }
    //printf("Latent %d: ", i);
    //for(j=0; j < 10; j++){
    //  printf("%d = %f, ", j, last.latent[j]);
    //}
    //printf("\n");
     
    // tau
    gibbs_tau(&last, &param, Y, spec, dims, 
              predictor, 
              K, N, n, nY, nfix, cumnfix, &totnfix, n_j, i_j, nk);
    if(whatsave[8]){
      for(j = 0; j < dimswithK[8]; j++){
        tau[i*dimswithK[8] + j] = last.tau[j];
      }
    }
    //printf("Tau %d: ", i);
    //for(j=0; j < dimswithK[8]; j++){
    //  printf("%d = %f, ", j, last.tau[j]);
    //}
    //printf("\n");
    
    // rantau
    gibbs_rantau(&last, &param, spec, dims,  
                 K, &totnY, nran, cumnran, &totnran);
    if(whatsave[10]){
      for(j = 0; j < dimswithK[10]; j++){
        rantau[i*dimswithK[10] + j] = last.rantau[j];
      }
    }
    //printf("RanTau %d: ", i);
    //for(j=0; j < dimswithK[10]; j++){
    //  printf("%d = %f, ", j, last.rantau[j]);
    //}
    //printf("\n");
    
    
    // beta 
    gibbs_beta(&last, &param, Id, Y, X, spec, dims, K, N, n, nY, 
               &totnY, nfix, nran, cumnfix, cumnran, &totnfix, &totnran, 
               FormulaF, FormulaR, n_j, i_j, nUk);
    if(whatsave[7]){
      for(j = 0; j < dimswithK[7]; j++){
        beta[i*dimswithK[7] + j] = last.beta[j];
      }
    }
    //printf("Beta %d: ", i);
    //for(j=0; j < dimswithK[7]; j++){
    //    printf("%f, ", last.beta[j]);
    //}
    //printf("\n");
    
    // mu
    gibbs_mu(&last, &param, spec, dims, K, n, &totnran, nUk);
    if(whatsave[9]){
      for(j = 0; j < dimswithK[9]; j++){
        mu[i*dimswithK[9] + j] = last.mu[j];
      }
    }
    //printf("Mu %d: ", i);
    //for(j=0; j < dimswithK[9]; j++){
    //    printf("%f, ", last.mu[j]);
    //}
    //printf("\n");
    
    // InvQ
    gibbs_InvQ(&last, &param, spec, calc, dims, dimswithK, 
               K, &totnran, &i, Q, detQ, detInvQ);
    if(whatsave[12]){
      for(j = 0; j < dimswithK[12]; j++){
        InvQ[i*dimswithK[12] + j] = last.InvQ[j];
      }
    }

    // InvSigma
    gibbs_InvSigma(&last, &param, spec, calc, dims, dimswithK, 
                   K, n, &totnran, &i, nUk, 
                   Sigma, detSigma, detInvSigma, sdSigma, corSigma);
    if(whatsave[11]){
      for(j = 0; j < dimswithK[11]; j++){
        InvSigma[i*dimswithK[11] + j] = last.InvSigma[j];
      }
    }
    //printf("InvSigma %d: ", i);
    //for(j=0; j < dimswithK[11]; j++){
    //  printf("%f, ", last.InvSigma[j]);
    //}
    //printf("\n");
    
    // b
    gibbs_b(&last, &param, Y, X, spec, dims, FormulaF, FormulaR, 
            N, n, K, ncat, nY, &totnY, nfix, nran, cumnfix, cumnran, 
            &totnfix, &totnran, &i, n_j, &max_n_j, i_j, ZtZ, dimZtZ);
    if(whatsave[13]){
      for(j = 0; j < dims[13]; j++){
        b[i*dims[13] + j] = last.b[j];
      }
    }
    //printf("Random effects %d: ", i);
    //for(j=0; j < 10; j++){
    //  printf("%f, ", last.b[j]);
    //}
    //printf("\n");

    // pUik_nonb... not in last. --> must be saved within this function 
    if(calc[0]){
      gibbs_pUik_nonb(&last, &param, Y, X, spec, dims, FormulaF, FormulaR, 
                      N, n, K, ncat, nY, &totnY, nfix, nran, cumnfix, cumnran, 
                      &totnfix, &totnran, &i, n_j, &max_n_j, i_j, isnonzero, logp_base, 
                      nUk, nk, pUik_nonb, ZtZ, dimZtZ);
      // results should be stored after last i iterations 
    }
  } // end of for i 
  
  printf("\n");
  fflush(stdout);
  
  /*
   * Transferring last state as an output
   */
  
  /*
  
  NOT necessary - pointers in the structures last are set to the memory space given by function
  
  // last U
  for(i = 0; i < dimswithK[1]; i++){
    last_U[i] = last.U[i];
  }
  
  // last other parameters
  int counter = 0;
  // w 
  for(i = 0; i < dimswithK[0]; i++){
    veclast[counter++] = last.w[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[2]; i++){
    veclast[counter++] = last.pUik[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[3]; i++){
    veclast[counter++] = last.gamma[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[6]; i++){
    veclast[counter++] = last.latent[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[7]; i++){
    veclast[counter++] = last.beta[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[8]; i++){
    veclast[counter++] = last.tau[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[9]; i++){
    veclast[counter++] = last.mu[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[10]; i++){
    veclast[counter++] = last.rantau[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[11]; i++){
    veclast[counter++] = last.InvSigma[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[12]; i++){
    veclast[counter++] = last.InvQ[i]; // counter is increased after copying the value
  }
  for(i = 0; i < dimswithK[13]; i++){
    veclast[counter++] = last.b[i]; // counter is increased after copying the value
  }
  
  */
  
}




