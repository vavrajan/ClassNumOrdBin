/* Predictor calculation */
#include <R.h>
#include <Rmath.h>
#include <stdio.h>

void calculate_predictor( 
    /** IN parameters **/
    int* Id,              // [N] IDs
    double* X,      // [N*#regr]    ID + regressors
    int* spec,     // [7]                class-specific parameters
    // order:   [0] gamma, 
    //          [1] tau, 
    //          [2] beta, 
    //          [3] rantau, 
    //          [4] mu, 
    //          [5] InvQ, 
    //          [6] InvSigma
    /** Output **/
    double* predictor,      //          [*N * totnY] current value of predictor
    //int* posun,
    //int* isspec, 
    /** Parameters describing dimensions **/
    int* FormulaF,    // [sum(nY)]  numbers of columns of X that should be used for FIXED  effects of modelled responses
    int* FormulaR,    // [sum(nY)]  numbers of columns of X that should be used for RANDOM effects of modelled responses
    int* N,           // [1]        total number of observations
    int* n,           // [1]        total number of subjects (different ids)
    int* dims,        // [23]       the length of subarray that corresponds to one state (disected by various parameters)
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
    int* totnY,       // [1]        total number of responses
    int* nfix,        // [sum(nY)]  number of FIXED  regressors for each response
    int* nran,        // [sum(nY)]  number of RANDOM regressors for each response
    int* cumnfix,     // [sum(nY)]  cummulative sums of nfix values and -1 (for C indexing)
    int* cumnran,     // [sum(nY)]  cummulative sums of nran values and -1 (for C indexing)
    /** Arrays with necessary parameters from one state **/
    double* beta,         // [totnfix(* K)]   Beta parameters for fixed effects
    double* b,            // [totnran * n]    Random effects
    int* U             // [n]              What is the current class of i-th subject
)
{
  /*** Declarations ***/
  /*** ------------ ***/
  int y;                // number of current response
  int i1, i2;           // looping indeces
  int cat;              // number of class (category)
  int id;               // number of ID
  int posbeta;          // position in beta parameter - depends on the current class
  int xdim;
  int bdim;
  
  /*** Predictor calculation ***/
  /*** --------------------- ***/
  
  
  // Take the number of observation
  for(i1 = 0; i1 < *N; i1++){
    
    // What is the class this observation belongs to?
    id = Id[ i1 ];            
    cat = U[ id ];  
    
    // Is beta parameter class-specific?
    if(spec[2]){
      posbeta = cat;        // to shift at place corresponding to cat(th) beta value
      //*isspec = 1;
    }else{
      posbeta = 0;              // there is no shift needed - 
      // there is only one value of beta parameter
      //*isspec = 0;
    }
    //posun[i1] = posbeta;
    
    // Compute for each response (Numeric, Ordinal and Binary)  
    for(y = 0; y < *totnY; y++){
      
      predictor[i1 + *N * y] = 0;
      
      // fixed effects
      for(i2 = 0; i2 < nfix[y]; i2++){
        xdim = i1 + *N * FormulaF[i2 + cumnfix[y]];
        bdim = i2 + cumnfix[y] + posbeta * dims[7];
        // add i2-th regressor for y multiplied by corresponding beta
        predictor[i1 + *N * y] += X[xdim] * beta[bdim];
      }
      
      // random effects
      for(i2 = 0; i2 < nran[y]; i2++){
        xdim = i1 + *N * FormulaR[i2 + cumnran[y]];
        bdim = id + *n * (cumnran[y] + i2);
        // add i2-th regressor for y multiplied by corresponding b
        predictor[i1 + *N * y] += X[xdim] * b[bdim];
      }
      
    } // end of for(y=0;...)
  } // end of for(i1 = 0)
  
} // end of function calculate_predictor
