/*
 * Generating latent variables for ordinal and binary outcomes
 */


#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>

#include "myrtruncnorm.h"

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
                )    
{ 
  /** Declarations **/
  int l, m;               // looping indeces
  int y;                  // variable index
  double bounds[*maxnbounds]; // array for bounds for gamma parameter (may be longer than actually needed)
  int dimgammay;          // dimension of gamma parameter for specific response y
  // calculated as #(categories) - 1 - 1(for counds and fixation of the first one) 
  int cat;                // category
  int cumdimgammay = 0;   // cummulated dimension of gamma
  double sd = 1.0;        // to be passed as adress to a function
  
  /* 
  * Generating latent variables and gammas for Ordinal responses 
  */
  
  for(y = 0; y < nY[1]; y++){
    // for cycle through ordinal variables
    // there is only one gamma parameter for each ordinal response
      
    dimgammay = ncat[y]-2;
    
    // setting bounds
    bounds[0] = gamma1[y];
    for(m = 1; m <= dimgammay; m++){
      bounds[m] = pgamma[m - 1 + cumdimgammay];
    }
    
    for(l = 0; l < n_j[*j]; l++){
      cat = Y[i_j[*j + *nnew * l] + *Nnew * (y + nY[0])]; // conversion of type double into type int
      // cat in {0, 1, 2, ..., ncat[y]-1}
      // generate new latent
      rtruncnorm(latent + l + n_j[*j] * y,
                 &cat,
                 bounds,
                 &dimgammay,
                 predictor + l + n_j[*j] * (y + nY[0]),
                 &sd);
    }
    
    // increase cumdimgammay for next ordinal response
    cumdimgammay += dimgammay;
  } // end of for y in Ord
  
  
  /* 
   * Generating latent variables for Binary responses
   */ 
  dimgammay = 0;
  
  for(y = 0; y < nY[2]; y++){
    for(l = 0; l < n_j[*j]; l++){
      cat = Y[i_j[*j + *nnew * l] + *Nnew * (y + nY[0] + nY[1])]; // conversion of type double into type int
      // cat is either 0 or 1
      // bound dividing these two categories is just one given by hyperparameter gammaBin
      
      // generate new latent
      rtruncnorm(latent + l + n_j[*j] * (y + nY[1]),
                 &cat,
                 gammaBin,
                 &dimgammay,
                 predictor + l + n_j[*j] * (y + nY[0] + nY[1]),
                 &sd);
    }
  }  
  
  
} // end of void 
