/*
 * Generating from Truncated normal distribution
 */


#include <R.h>
#include <Rmath.h>

void rtruncnorm(double* gen,    // OUT [1] generated values
                int* cat,       // IN  [1] number of category
                double* bounds, // IN  [dimgammay+1] bounds:
                                // cat=0 | bounds[0] | cat=1 | bounds[1] | ... | bounds[dimgamma] | cat=dimgamma+1
                int* dimgammay, // IN  [1] how many bounds are randomized (1 is always fix)... dimgammay + 2 = total number of categories
                double* mean,   // IN  [1] mean value
                double* sd      // IN  [1] standard deviation
                )    
{ 
  /** Declarations **/
  double phia, phib;  // sum of all gamma distributed variables
  double x;           // generated value
  double u;           // generated uniformely distributed value
  
  if(*cat == 0){
    // we generate from N(mean,sd)|(-Inf, bounds[0]) 
    phia = 0.0;
    phib = pnorm((bounds[0]-*mean)/(*sd), 0.0, 1.0, 1, 0);

  }else if(*cat == *dimgammay+1){
    // we generate from N(mean,sd)|(bounds[*dimgamma], Inf)
    phia = pnorm((bounds[*dimgammay]-*mean)/(*sd), 0.0, 1.0, 1, 0);
    phib = 1.0;
    
  }else{
    // we generate from N(mean,sd)|(bounds[*cat-1], bounds[*cat])
    phia = pnorm((bounds[*cat-1]-*mean)/(*sd), 0.0, 1.0, 1, 0);
    phib = pnorm((bounds[*cat]-*mean)/(*sd),   0.0, 1.0, 1, 0);
    
  }
  
  u = runif(0.0, 1.0);
  x = phia + u * (phib - phia);
  x = qnorm(x, 0.0, 1.0, 1, 0);
  *gen = *sd * x + *mean;
  
}

