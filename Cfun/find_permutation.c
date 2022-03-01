


#include <R.h>
#include <Rmath.h>
#include <stdio.h>


/*** ===================================================================== ***/
/*** While cycle for finding the optimal permutations ***/
/*** ===================================================================== ***/

void FindPermutations(
      int* K,         // IN  [1]     number of clusters
      int* M,         // IN  [1]     number of iterations
      int* n,         // IN  [1]     number of subjects
      int* nperms,    // IN  [1]     K! - number of permutations
      int* perms,     // IN  [K!*K]  rows contain all permutations 
      double* probs,  // IN  [M*K*n] generated probabilities
      int* iteration, // OUT [1]     number of iterations needed to converge
      int* max_iter,  // IN  [1]     maximal number of iterations
      int* nu,        // OUT [M*K]   m-th row contains permutation of m-th iteration
      double* Q,      // OUT [K*n]   i-th column contains Q probabilities of i-th subject
      double* eps     // IN  [1]     tolerance for convergence
)
{
  /* Declarations */ 
  int converged = 0;
  //double newnu[*M * *K];
  //double newQ[*K * *n];
  double sum;
  double min;
  double Fnorm_Q;
  double Fnorm_nu;
  
  // indices
  int i, j, k, t;
  long coord;
  long pcoord;
  long nucoord;
  int theperm;
  
  printf("Before while cycle. \n");
  fflush(stdout);
  
  while(!converged){
    (*iteration)++;
    
    printf("Starting iteration %d \n", *iteration);
    fflush(stdout);
    // Step 1 - computing Q
    
    coord = 0;
    //pcoord = 0;
    Fnorm_Q = 0.0;
    for(i = 0; i <*n; i++){
      nucoord = 0;
      for(k = 0; k <*K; k++){
        //newQ[coord] = 0.0;
        sum = 0.0;
        for(t = 0; t <*M; t++){
          // nucoord = t + *M * k --> index of t-th permutation of k-th element within matrix nu
          //printf("Step 1: i=%d, k=%d, t=%d, coord=%ld, nucoord=%ld, pcoord=%ld, sum=%f\n", 
          //       i, k, t, coord, nucoord, pcoord, sum);
          pcoord = t + nu[nucoord] * *M + i * *K * *M;
          //printf("Step 1: i=%d, k=%d, t=%d, coord=%ld, nucoord=%ld, pcoord=%ld, sum=%f\n", 
          //       i, k, t, coord, nucoord, pcoord, sum);
          //newQ[coord] += probs[pcoord];
          sum += probs[pcoord];
          //printf("Step 1: i=%d, k=%d, t=%d, coord=%ld, nucoord=%ld, pcoord=%ld, sum=%f\n", 
          //       i, k, t, coord, nucoord, pcoord, sum);
          //fflush(stdout);
          nucoord++; // another t in given k
          
        }
        //newQ[coord] /= *M; // making mean out of a sum
        sum /= *M; // making mean out of a sum
        //Fnorm_Q += (Q[coord] - newQ[coord])*(Q[coord] - newQ[coord]);
        Fnorm_Q += (Q[coord] - sum)*(Q[coord] - sum);
        //printf("Step 1: i=%d, k=%d, coord=%ld, sum=%f, Fnorm_Q=%f\n", 
        //              i, k, coord, sum, Fnorm_Q);
        Q[coord] = sum;
        coord++; // another k in given i
      }
    }
    
    printf("End of Step 1: iteration = %d, Fnorm_Q = %f\n", *iteration, Fnorm_Q);
    fflush(stdout);
    
    // Step 2 - optimizing over permutations
    Fnorm_nu = 0.0;
    for(t = 0; t < *M; t++){
      theperm = 0;
      min = 100000000000.0;
      for(j = 0; j <*nperms; j++){
        sum = 0.0;
        coord = 0;
        for(i = 0; i <*n; i++){
          for(k = 0; k <*K; k++){
            pcoord = t + perms[j + *nperms * k] * *M + i * *K * *M;
            //printf("Step 2: t=%d, j=%d, i=%d, k=%d, coord=%ld, pcoord=%ld, min=%f, theperm=%d, Q=%f\n", 
            //       t, j, i, k, coord, pcoord, min, theperm, Q[coord]);
            sum += probs[pcoord] * (log(probs[pcoord]) - log(Q[coord]));
            //printf("Step 2: t=%d, j=%d, i=%d, k=%d, coord=%ld, pcoord=%ld, sum=%f\n", 
            //       t, j, i, k, coord, pcoord, sum);
            //fflush(stdout);
            coord++;
          }
        }
        if(min > sum){
          min = sum;
          theperm = j;
        }
      }
      
      for(k = 0; k < *K; k++){
        //newnu[t + *M * k] = perms[theperm + *nperms * k];
        coord = theperm + *nperms * k;
        nucoord = t + *M * k;
        Fnorm_nu += (nu[nucoord] - perms[coord])*(nu[nucoord] - perms[coord]);
        nu[nucoord] = perms[coord];
      }
    }
    printf("End of Step 2: iteration = %d, Fnorm_nu = %f\n", *iteration, Fnorm_nu);
    fflush(stdout);
    
    // Did it converge?
    // I define convergence as Frobenius norm of (Q-newQ) < epsilon
    //Fnorm_Q = 0.0;
    //for(j = 0; j < *n * *K; j++){
    //  Fnorm_Q += (Q[j] - newQ[j])*(Q[j] - newQ[j]);
    //}
    Fnorm_Q = sqrt(Fnorm_Q);
    Fnorm_nu = sqrt(Fnorm_nu);
    converged = ((Fnorm_Q < *eps) || (Fnorm_nu < *eps) || (*iteration == *max_iter));
    printf("iteration = %d, Fnorm_Q = %f, Fnorm_nu = %f\n", *iteration, Fnorm_Q, Fnorm_nu);
    fflush(stdout);
    
    // Saving last values
    //for(j = 0; j < *M * *K; j++){
    //  nu[j] = newnu[j];
    //}
    //for(j = 0; j < *n * *K; j++){
    //  Q[j] = newQ[j];
    //}
    
  }
}