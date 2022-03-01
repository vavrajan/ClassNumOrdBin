#ifndef GIBBS_PUIK_H
#define GIBBS_PUIK_H

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
);


#endif