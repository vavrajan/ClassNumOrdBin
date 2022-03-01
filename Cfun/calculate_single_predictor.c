/* Predictor calculation for single subject */
#include <R.h>
#include <Rmath.h>
#include <stdio.h>

void calculate_single_predictor( 
    /** IN parameters **/
    double* X,            // [Nnew*#regr]    ID + regressors
    int* Nnew,            // [1] total number of observations
    int* nnew,            // [1] total number of subjects
    /** Output **/
    double* predictor,    // [*N * totnY] current value of predictor
    /** Subject **/
    int* j,               // [1] subject number
    int* n_j,             // [*nnew]   number of observations dedicated to j-th subject
    int* i_j,             // [*nnew * max_n_j] matrix of indeces to j-th subject (j-th row = j-th subject)
    /** Parameters describing dimensions **/
    int* FormulaF,        // [sum(nY)]  numbers of columns of X that should be used for FIXED  effects of modelled responses
    int* FormulaR,        // [sum(nY)]  numbers of columns of X that should be used for RANDOM effects of modelled responses
    int* totnY,           // [1]        total number of responses
    int* nfix,            // [sum(nY)]  number of FIXED  regressors for each response
    int* nran,            // [sum(nY)]  number of RANDOM regressors for each response
    int* cumnfix,         // [sum(nY)]  cummulative sums of nfix values and -1 (for C indexing)
    int* cumnran,         // [sum(nY)]  cummulative sums of nran values and -1 (for C indexing)
    /** Arrays with necessary parameters from one state **/
    double* beta,         // [totnfix]   Beta parameters for fixed effects
    double* b             // [totnran * n]    Random effects
)
  
{   
  // Declarations
  int y;                // number of current response
  int pos;              // position in predictor array
  int i1, i2;           // looping indeces
  int xdim;             // position in X matrix
  int bdim;             // position in b or beta
  
  for(i1 = 0; i1 < n_j[*j]; i1++){
    // for all observations dedicated to j-th subject
    
    // Compute for each response (Numeric, Ordinal and Binary)  
    for(y = 0; y < *totnY; y++){
      
      pos = i1 + n_j[*j] * y; // row i1, column y
      predictor[pos] = 0;
      
      // fixed effects
      for(i2 = 0; i2 < nfix[y]; i2++){
        bdim = i2 + cumnfix[y];
        xdim = i_j[*j + *nnew * i1]  + *Nnew * FormulaF[bdim];
        // add i2-th regressor for y multiplied by corresponding beta
        predictor[pos] += X[xdim] * beta[bdim];
      }
      
      // random effects
      for(i2 = 0; i2 < nran[y]; i2++){
        bdim = i2 + cumnran[y];
        xdim = i_j[*j + *nnew * i1]  + *Nnew * FormulaR[bdim];
        // add i2-th regressor for y multiplied by corresponding beta
        predictor[pos] += X[xdim] * b[bdim];
      }
      
    } // end of for(y=0;...)
  }
  
}