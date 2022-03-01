#ifndef CALCULATE_SINGLE_PREDICTOR_H
#define CALCULATE_SINGLE_PREDICTOR_H

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
    );

#endif