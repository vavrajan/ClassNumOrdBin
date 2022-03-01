#ifndef GIBBS_LATENT_H
#define GIBBS_LATENT_H

void gibbs_latent(struct str_state* last,  // OUT+IN last known values of generated parameters
                  struct str_param* param, // IN hyperparameters
                  double* Y,            // IN [*N*totnY]
                  double* predictor,    // IN [*N * totnY] current value of predictor
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
                  int* N,           // IN [1] number of observations 
                  int* n,           // IN [1] number of subjects
                  int* K,           // IN [1] number of classes
                  int* ncat,        // IN [nY[1 .. Ords]]  the counts of categories of ordinal variables
                  int maxnbounds,   // IN [1] maximal number of thresholds gamma for one variable
                  int* nY,          // IN [3]        counts of Nums, Ords and Bins variables
                  int* n_j,         // IN [n] number of observations dedicated to j-th subject
                  int* i_j,         // IN [n * max_n_j] indeces dedicated to j-th subject 
                  int* savemax,     // IN [1] should it save max?
                  int* savemin,     // IN [1] should it save min?
                  double* pmax,     // IN [1] pointer to saving max
                  double* pmin      // IN [1] pointer to saving min
);  

#endif