#ifndef TRUNC_NORM_H
#define TRUNC_NORM_H

void rtruncnorm(double* gen,    // OUT [1] generated values
                int* cat,       // IN  [1] number of category
                double* bounds, // IN  [dimgammay+1] bounds:
                // cat=0 | bounds[0] | cat=1 | bounds[1] | ... | bounds[dimgamma] | cat=dimgamma+1
                int* dimgammay, // IN  [1] how many bounds are randomized (1 is always fix)... dimgammay + 2 = total number of categories
                double* mean,   // IN  [1] mean value
                double* sd      // IN  [1] standard deviation
);   


#endif