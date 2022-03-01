#ifndef GIBBS_INVSIGMA_H
#define GIBBS_INVSIGMA_H

void gibbs_InvSigma(struct str_state* last,  // OUT+IN last known values of generated parameters
                    struct str_param* param, // IN hyperparameters
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
                    int* dimswithK,       // [23]       the length of subarray that corresponds to one state 
                    //            (disected by various parameters, also multiplication by K incorporated when such parameters is class-specific)
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
                    int* i,                 // IN [1] iteration number
                    int* nUk,               // IN [K] number of subjects     in classes 1, ..., K
                    double* Sigma,          // [(B+M)*totnran*(totnran+1)/2(* K)] Prior variance matrix of random effects
                    double* detSigma,       // [(B+M)(* K)] determinant of Sigma
                    double* detInvSigma,    // [(B+M)(* K)] determinant of InvSigma
                    double* corSigma,       // [(B+M)*totnran*(totnran-1)/2(* K)] Prior correlation matrix of random effects
                    double* sdSigma         // [(B+M)*totnran(* K)] Prior standard deviations of random effects
);  

#endif