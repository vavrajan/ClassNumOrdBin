#ifndef GIBBS_MU_H
#define GIBBS_MU_H

void gibbs_mu(struct str_state* last,  // OUT+IN last known values of generated parameters
              struct str_param* param, // IN hyperparameters
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
              int* K,                 // IN [1] number of classes
              int* n,                 // IN [1] number of subjects
              int* totnran,           // IN [1] total number of random regressors
              int* nUk                // IN [K] number of subjects     in classes 1, ..., K
); 

#endif