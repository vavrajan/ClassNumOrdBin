#############################################
###   Compile C functions into ".dll"     ###
###---------------------------------------###

### Purpose of this script is to compile used C functions
### into ".dll" file that is loaded within R functions

### First set up the working directory
### ROOT should contain the path to directory where (ending with "/")
### the following directories have been saved to
### - ./Rfun - contains all R functions
### - ./Cfun - contains all C functions

ROOT <- "C:/Users/Vavra/Desktop/SKOLA/PhD studium/Disertacka/Classification_final/"
### set the paths to Cfun and Rfun directories
CROOT <- paste0(ROOT, "Cfun/")
RROOT <- paste0(ROOT, "Rfun/")


### Compiling C functions
### libraries "withr" and "callr" needed
library("withr")
library("callr")

setwd(CROOT)

# Compile all functions needed for MCMC sampling
# creates "class_num_ord_bin.dll" that is used within R function "GibbsClassNumOrdBin"
# make sure that the module is not loaded elsewhere
out <- rcmd(cmd = "SHLIB", cmdargs = c("class_num_ord_bin.c",
                                       "calculate_predictor.c",
                                       "cholesky.c",
                                       "matrmult.c",
                                       "myrdirichlet.c",
                                       "myrtruncnorm.c",
                                       "myrwishart.c",
                                       "gibbs_b.c",
                                       "gibbs_beta.c",
                                       "gibbs_InvQ.c",
                                       "gibbs_InvSigma.c",
                                       "gibbs_latent.c",
                                       "gibbs_mu.c",
                                       "gibbs_pUik.c",
                                       "gibbs_pUik_nonb.c",
                                       "gibbs_rantau.c",
                                       "gibbs_tau.c",
                                       "structures.h"))
cat(out$stderr)
out



# Compile all functions needed for calculation of classification probabilities
# creates "pUnewk_calculation.dll" that is used within R function "pUnewk_pmvnorm" 
# make sure that the module is not loaded elsewhere
out <- rcmd(cmd = "SHLIB", cmdargs = c("pUnewk_pmvnorm.c",
                                       "calculate_single_predictor.c",
                                       "cholesky.c",
                                       "matrmult.c",
                                       "myrtruncnorm.c",
                                       "rmvnorm.c",
                                       "gen_latent.c"
))
cat(out$stderr)
out

# single "find_permutation" function which performs procedure by Stevens (2000) 
# for detection of label-switching problem and suggestion of its remedy
out <- rcmd(cmd = "SHLIB", cmdargs = "find_permutation.c")
cat(out$stderr)
out
