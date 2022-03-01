/*
 * Generating gamma and latent variables for ordinal and binary outcomes
 */


#include <R.h>
#include <Rmath.h>
#include <math.h>
#include <stdio.h>

#include "myrtruncnorm.h"

#include "structures.h"


#define LARGE_VALUE 1000000000.0
#define LOW_VALUE  -1000000000.0

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
                  )    
{ 
  /** Declarations **/
  int j, k, l, m;         // looping indeces
  int y;                  // variable index
  double bounds[maxnbounds]; // array for bounds for gamma parameter (may be longer than actually needed)
  int dimgammay;          // dimension of gamma parameter for specific response y
                          // calculated as #(categories) - 1 - 1(for counds and fixation of the first one) 
  int cat;                // category
  double maxs[maxnbounds];   // array of max values (may be longer than actually needed)
  double mins[maxnbounds];   // array of min values (may be longer than actually needed)
  int cumdimgammay = 0;   // cummulated dimension of gamma
  double sd = 1.0;        // to be passed as adress to a function
  
  /* 
   * Generating latent variables and gammas for Ordinal responses 
   */
  
  for(y = 0; y < nY[1]; y++){
    // for cycle through ordinal variables
    // Is gamma parameter class-specific?
    if(spec[0]){
      // we need to generate new K versions of gamma parameter
      for(k = 0; k < *K; k++){
        dimgammay = ncat[y]-2;
        
        // setting bounds
        bounds[0] = (*param).gamma1[y];
        for(m = 1; m <= dimgammay; m++){
          bounds[m] = (*last).gamma[m - 1 + cumdimgammay + k * dims[3]];
        }
        
        // setting maxs and mins
        for(m = 0; m < dimgammay; m++){
          maxs[m] = LOW_VALUE;
          mins[m] = LARGE_VALUE;
        }
        
        // generating new latent variables
        for(j = 0; j < *n; j++){
          // j is the subject
          if((*last).U[j] == k){
            // generate only if j-th cubject lies in class k
            for(l = 0; l < n_j[j]; l++){
              // l-th observation of subject j
              cat = Y[i_j[j + *n * l] + *N * (y + nY[0])]; // conversion of type double into type int
              // generate new latent
              rtruncnorm((*last).latent + i_j[j + *n * l] + *N * y,
                         &cat,
                         bounds,
                         &dimgammay,
                         predictor + i_j[j + *n * l] + *N * (y + nY[0]),
                         &sd);
              // updating maxs
              if(cat > 0 && cat <= dimgammay && maxs[cat-1] < (*last).latent[i_j[j + *n * l] + *N * y]){
                maxs[cat-1] = (*last).latent[i_j[j + *n * l] + *N * y];
              }
              // updating mins
              if(cat > 1 && cat <= dimgammay+1 && mins[cat-2] > (*last).latent[i_j[j + *n * l] + *N * y]){
                mins[cat-2] = (*last).latent[i_j[j + *n * l] + *N * y];
              }
            }
          }
        } // end of for j in subjects
      
      // save max values  
      if(*savemax){
        for(m = 0; m < dimgammay; m++){
          pmax[m + cumdimgammay + k*dims[5]] = maxs[m];
        }
      }

      // save min values
      if(*savemin){
        for(m = 0; m < dimgammay; m++){
          pmin[m + cumdimgammay + k*dims[4]] = mins[m];
        }
      }
      
      // generate new gamma
      for(m = 0; m < dimgammay; m++){
        (*last).gamma[m + cumdimgammay + k*dims[3]] = runif(maxs[m], mins[m]);
      }
        
    } // end of for k in classes
      
    }else{
      // there is only one gamma parameter for each ordinal response

      dimgammay = ncat[y]-2;
      
      // setting bounds
      bounds[0] = (*param).gamma1[y];
      for(m = 1; m <= dimgammay; m++){
        bounds[m] = (*last).gamma[m - 1 + cumdimgammay];
      }
      
      /*  printf("Bounds: ");
        for(m = 0; m<= dimgammay; m++){
        printf("%f, ", bounds[m]);
        }
        printf("\n");
       */
      
      
      // setting maxs and mins
      for(m = 0; m < dimgammay; m++){
        maxs[m] = LOW_VALUE;
        mins[m] = LARGE_VALUE;
      }
      
      // generating new latent variables
      for(j = 0; j < *N; j++){
        // j is the j-th observation
        
        cat = Y[j + *N * (y + nY[0])]; // conversion of type double into type int
        // cat in {0, 1, 2, ..., ncat[y]-1}
        // generate new latent
        rtruncnorm((*last).latent + j + *N * y,
                   &cat,
                   bounds,
                   &dimgammay,
                   predictor + j + *N * (y + nY[0]),
                   &sd);
        // updating maxs
        //if(j < 10){
        //  printf("maxTF = %d, minTF = %d\n",
        //      (cat > 0) && (cat <= dimgammay) && (maxs[cat-1] < (*last).latent[j + *N * y]),
        //      (cat > 1) && (cat <= (dimgammay+1)) && (mins[cat-2] > (*last).latent[j + *N * y]));
        //}
        if((cat > 0) && (cat <= dimgammay) && (maxs[cat-1] < (*last).latent[j + *N * y])){
          maxs[cat-1] = (*last).latent[j + *N * y];
        }
        // updating mins
        if((cat > 1) && (cat <= (dimgammay+1)) && (mins[cat-2] > (*last).latent[j + *N * y])){
          mins[cat-2] = (*last).latent[j + *N * y];
        }
        //if(j < 10){
        //  printf("cat = %d, cat-2 = %d, latent = %f, maxTF = %d, minTF = %d\n", 
        //   cat, cat-2, (*last).latent[j + *N * y], 
        //   (cat > 0) && (cat <= dimgammay) && (maxs[cat-1] < (*last).latent[j + *N * y]),
        //   (cat > 1) && (cat <= (dimgammay+1)) && (mins[cat-2] > (*last).latent[j + *N * y]));
        //  printf("Maxs = ");
        //  for(m = 0; m<dimgammay;m++){
        //    printf("%f", maxs[m]);
        //  }
        //  printf("\nMins = ");
        //  for(m = 0; m<dimgammay;m++){
        //    printf("%f", mins[m]);
        //  }
        //  printf("\n");
        //  }
        
      } // end of for j in observations
      
      // save max values  
      if(*savemax){
        for(m = 0; m < dimgammay; m++){
          pmax[m + cumdimgammay] = maxs[m];
        }
      }
      
      // save min values
      if(*savemin){
        for(m = 0; m < dimgammay; m++){
          pmin[m + cumdimgammay] = mins[m];
        }
      }
      
      // generate new gamma
      for(m = 0; m < dimgammay; m++){
        (*last).gamma[m + cumdimgammay] = runif(maxs[m], mins[m]);
      }
    } // end of else of if spec["gamma"]
    
    // increase cumdimgammay for next ordinal response
    cumdimgammay += dimgammay;
  } // end of for y in Ord
  
  
  /* 
   * Generating latent variables for Binary responses
   */ 
  dimgammay = 0;
    
  for(y = 0; y < nY[2]; y++){
    for(j = 0; j < *N; j++){
      cat = Y[j + *N * (y + nY[0] + nY[1])]; // conversion of type double into type int
      // cat is either 0 or 1
      // bound dividing these two categories is just one given by hyperparameter gammaBin
      
      // generate new latent
      rtruncnorm((*last).latent + j + *N * (y + nY[1]),
                 &cat,
                 (*param).gammaBin,
                 &dimgammay,
                 predictor + j + *N * (y + nY[0] + nY[1]),
                 &sd);
    }
  }  

  
} // end of void 
