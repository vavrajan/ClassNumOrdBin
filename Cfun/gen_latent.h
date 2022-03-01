#ifndef GEN_LATENT_H
#define GEN_LATENT_H

void gen_latent(double* latent,     // OUT [*Nnew * (nY[1]+nY[2])]
                double* predictor,  // IN [n_j[j] * totnY] predictor values for j-th subject
                double* Y,          // IN [*N*totnY]
                int* Nnew,          // IN [1] number of observations 
                int* nnew,          // IN [1] number of subjects
                int* ncat,          // IN [nY[1 .. Ords]]  the counts of categories of ordinal variables
                int* nY,            // IN [3]        counts of Nums, Ords and Bins variables
                int* j,             // IN [1] subject number
                int* n_j,           // IN [n] number of observations dedicated to j-th subject
                int* i_j,           // IN [n * max_n_j] indeces dedicated to j-th subject
                int* maxnbounds,    // IN [1] maximal number of thresholds gamma for one variable
                double* pgamma,     // IN [dims[2]] gamma parameter
                double* gamma1,     // [nY[1]]  fixed thresholds for     Ordinal variables
                double* gammaBin    // [1]      fixed threshold  for all Binary  variables
                );
#endif