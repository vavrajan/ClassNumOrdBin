#ifndef GIBBS_INVQ_H
#define GIBBS_INVQ_H

void gibbs_InvQ(struct str_state* last,  // OUT+IN last known values of generated parameters
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
                
                int* K,                 // IN  [1] number of classes
                int* totnran,           // IN [1] total number of random regressors
                int* i,                 // IN [1] iteration number
                double* Q,              // [(B+M)*totnran*(totnran+1)/2(* K)] Prior scale matrix of InvSigma
                double* detQ,           // [(B+M)(* K)] determinant of Q
                double* detInvQ         // [(B+M)(* K)] determinant of InvQ
); 

#endif