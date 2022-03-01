################################################
###   Tutorial use of GibbsClassNumOrdBin    ###
###------------------------------------------###

### Purpose of this script is to demonstrate the use of implemented functions
### for estimating model for longitudinal data of numeric, ordinal and binary nature
### for details of methodology, first read 
### Classification Based on Multivariate Mixed Type Longitudinal Data 
### with an application to the EU-SILC database
### by Jan Vávra and Arnošt Komárek in ADAC 2022

### First set up the working directory
### ROOT should contain the path to directory where (ending with "/")
### the following directories have been saved to
### - ./Rfun - contains all R functions
### - ./Cfun - contains all C functions

ROOT <- "path/to/your/directory/ClassNumOrdBin/"
### set the paths to subdirectories
CROOT <- paste0(ROOT, "Cfun/")  # C functions
RROOT <- paste0(ROOT, "Rfun/")  # R functions
RDATA <- paste0(ROOT, "RData/") # for saved RData files
FIG <- paste0(ROOT, "Figures/") # for saved figures
TAB <- paste0(ROOT, "Tables/")  # for saved tables

setwd(ROOT)

###----------------------------------------------------
### Loading needed libraries, R and C functions
###----------------------------------------------------

### needed R libraries
library("coda")
library("mvtnorm")
library("HDInterval")
library("colorspace")

### loading R functions
source(paste0(RROOT, "GenerateSimData.R"))
source(paste0(RROOT, "plotting_functions.R"))
source(paste0(RROOT, "function_GibbsClassNumOrdBin.R"))
source(paste0(RROOT, "FromCtoList.R"))
source(paste0(RROOT, "FromListtoC.R"))
source(paste0(RROOT, "FromCtoMatrix.R"))
source(paste0(RROOT, "FromMatrixtoC.R"))
source(paste0(RROOT, "ToCodaMCMC.R"))
source(paste0(RROOT, "permutation_functions.R"))
source(paste0(RROOT, "EachIdGradually.R"))
source(paste0(RROOT, "pUnewk_C_pmvnorm.R"))

### loading C functions (in Windows ".dll")
# if needed, you can recompile it using "compile_Cfun.R"
dyn.load(paste0(CROOT, "class_num_ord_bin.dll"))
dyn.load(paste0(CROOT, "pUnewk_pmvnorm.dll"))
dyn.load(paste0(CROOT, "find_permutation.dll"))


### loading C functions (".so" version)
#dyn.load(paste0(CROOT, "class_num_ord_bin.so"))
#dyn.load(paste0(CROOT, "pUnewk_pmvnorm.so"))
#dyn.load(paste0(CROOT, "find_permutation.so"))



###-------------------------------------------------------
### Generating sample longitudinal data
###-------------------------------------------------------

### loading true parameter values
load(paste0(RDATA, "sim_true_values.RData"))

n <- 1000           # sample size (number of subjects)
n_i <- 4            # observations per subject
K <- 3              # number of underlying clusters (works with K = 2 as well)
d <- "intercept"    # difference among cluster is in intercepts
# other possibilities are in
differences
r <- "intercept"    # random effects are formed by intercepts only
# other possibilities are in
randomparts
probs               # list of marginal probabilities of cluster allocation when generating
# we use "rep(1/K, K)" in K-th list

Betamu        # listed values for the true beta parameters in the following order
# [[K]][[k]][[r]][[d]][[Y]]
# [1] effect of Bernoulli(0.5) regressor
# [2] effect of Unif(0, 1) regressor
# [3] intercept (or mean of random intercept)
# [4] slope for time (or mean of slope random effect)

Eb            # listed values for mean of random effects
# [[K]][[k]][[r]][[d]]

Sigma         # listed values for covariance matrix of random effects
# [[K]][[k]][[r]][[d]]

tauNum        # precision of the numeric outcome
gam           # rest of gamma parameter for ordinal outcomes 
              # (first is -1 by default in GenerateData and later set in param below)
gamBin        # threshold for binary outcome
nY            # total number of outcomes
Nums          # names of Numeric outcomes
Ords          # names of Ordinal outcomes
Bins          # names of Binary outcomes


data <- GenerateData(n = n, n_i = n_i, 
                     K = K, difference = d, randompart = r,
                     probs = probs, Betamu = Betamu, Sigma = Sigma, Eb = Eb,
                     tauNum = tauNum, gam = gam, gamBin = gamBin,
                     nY = nY, Ords = Ords, Bins = Bins)

head(data, 12)

# plotting panel data using function PlotPanelDataIntoOne from "plotting_functions.R"
# fitted lines are from simple linear regression vs time
par(mar = c(4,4,1,1))
PlotPanelDataIntoOne(data, "Y1", K = K, subjects = 1:n)
PlotPanelDataIntoOne(data, "Y2", K = K, subjects = 1:n)
PlotPanelDataIntoOne(data, "fY2", K = K, addjitter = T, subjects = 1:n)
PlotPanelDataIntoOne(data, "Y3", K = K, subjects = 1:n)
PlotPanelDataIntoOne(data, "fY3", K = K, addjitter = T, subjects = 1:n)


###----------------------------------------------------
### Preparing input parameters for MCMC sampling
###----------------------------------------------------

### data preparation for Gibbs sampler
data <- transform(data, intercept = 1)

# outcome dataset
Y <- data[, c("Y1", "fY2", "fY3", "subject")]
summary(Y[[Nums[1]]]) # numeric outcome
summary(Y[[Ords[1]]]) # ordinal outcome as a factor with levels 0, 1, 2, .., L
summary(Y[["fY3"]]) # fY3 is a factor with levels 0, 1
# However function GibbsClassNumOrdBin expects binary as numeric 0, 1 values
Y$Y3 <- as.numeric(as.character(Y$fY3))
summary(Y[[Bins[1]]])

summary(Y[,Ys<-c(Nums, Ords, Bins)])

# covariates dataset
X <- data[, c("X1", "X2", "intercept", "time", "subject")]

### Auxiliary parameters for R functions independent of nsim
Id <- "subject" # name of the variable containing ids 

## list of formulas for fixed and random part of the model
## within the simulation all outcomes have the same formulas
## however, it can be changed to arbitrary 
## formula is given by a set of column names in X
## corresponding betas will then be in the order of these given names
Formula <- list()
# Fixed part
Formula[[Nums[1]]]$fixed <- Formula[[Ords[1]]]$fixed <- 
  Formula[[Bins[1]]]$fixed  <- switch(r,
                                      intercept = c("X1", "X2", "time"),
                                      slope = c("X1", "X2", "intercept"),
                                      both = c("X1", "X2"))
# if r=intercept, then only "X1", "X2" and "time" are fixed effects

# Random part
Formula[[Nums[1]]]$random <- Formula[[Ords[1]]]$random <- 
  Formula[[Bins[1]]]$random  <- switch(r,
                                       intercept = c("intercept"),
                                       slope = c("time"),
                                       both = c("intercept", "time"))
# if r=intercept, then only intercept forms the random effects 

## Specification of which parameters of the model should be latent group-specific
# follow the order of T/F values by the names given below
spec <- c(F, F, T, F, T, F, F)
names(spec) <- c("gamma", "tau", "beta", "rantau", "mu", "InvQ", "InvSigma")
# gamma = thresholding parameters for ordinal outcomes (RECOMMENDED TO USE ALWAYS FALSE, not tested with T)
# tau = precision parameter for the error terms of numeric outcomes
# beta = fixed effects
# rantau = prior precisions for mu (auxiliarty parameter)
# mu = mean of random effects (not forced to be centred!!!)
# InvQ = inverse of Q matrix used in Wishart prior for InvSigma
# InvSigma = inverse of matrix Sigma - covariance matrix of random effects

# for now beta and mu are the only cluster-specific parameters
# that is cluster can differ only in parameters describing the effects of covariates


## Other calculable parameters 
# selected (T) parameters will be calculated and saved in the sampled chains
calc <- c(F,
          F,F,F,
          T,F,F,T,T)
names(calc) <- c("pUik_nonb", 
                 "Q", "detQ", "detInvQ",
                 "Sigma", "detSigma", "detInvSigma", "sdSigma", "corSigma")
# pUik_nonb = P(U_i = k | all but b_i), random effects integrated out
#              ... not recommended to be used (just a remainder from development) !
#             use F always
# Q = original matrix Q, that is, inverse of InvQ
# detQ, detInvQ = determinant of respective matrices
# Sigma = covariance matrix of random effects, inverse of InvSigma
# detSigma, detInvSigma = determinants of respective matrices
# sdSigma = sqrt(diag(Sigma)) = standard deviations of random effects 
# corSigma = diag(sdSigma^{-1}) %*% Sigma %*% diag(sdSigma^{-1})
#          = correlation matrix of random effects



## What parameters should be saved?
# model works with a lot of parameters, some of which can be redundant to be saved
whatsave <- c(T,F,F,
              T,F,F,F,
              T,T,
              T,F,
              T,F,F)
names(whatsave) <- c("w", "U", "pUik",
                     "gamma", "min", "max", "latent",
                     "beta", "tau", 
                     "mu", "rantau",
                     "InvSigma", "InvQ", "b")
# w = P(U_i = k) = marginal cluster allocation probabilities
# U = current allocation of each of n subjects
# pUik = full-conditional probabilities used to sample new U within Gibbs procedure
# gamma = thresholds for ordinal outcomes
# min, max = interval bound used when sampling new gamma parameter
# latent = latent numeric outcomes for categorical outcomes
# beta = fixed effects
# tau = precisions of numeric outcomes
# mu = mean of random effects
# rantau = prior precision of random effects
# InvSigma = inverse of Sigma matrix
# InvQ = inverse of Q matrix (matrix used in Wishart prior for Sigma)
# b = random effects of all subjects

### Dimensions
nY <- sapply(list(Nums, Ords, Bins), length)
names(nY) <- c("Nums", "Ords", "Bins")
nfix <- nran <- numeric(sum(nY))
names(nfix) <- names(nran) <- c(Nums, Ords, Bins)
fixnames <- rannames <- numeric()
lfixnames <- lrannames <- list()
for(y in c(Nums, Ords, Bins)){
  nfix[y] <- length(Formula[[y]]$fixed)
  nran[y] <- length(Formula[[y]]$random)
  lfixnames[[y]] <- paste(y, Formula[[y]]$fixed, sep = "_")
  lrannames[[y]] <- paste(y, Formula[[y]]$random, sep = "_")
  fixnames <- c(fixnames, lfixnames[[y]])
  rannames <- c(rannames, lrannames[[y]])
}
nfix # d^F_r dimension of fixed effects separately by outcomes
lfixnames # list of names for the fixed effects
totnfix <- sum(nfix) # d^F dimension of fixed effects 

nran # d^F_r dimension of random effects separately by outcomes
lrannames # list of names for the random effects
totnran <- sum(nran) # d^R dimension of random effects


### Hyperparameters setting
param <- list(nu_0 = totnran + 1, # Wishart prior for InvSigma
              nu_1 = totnran + 1, # Wishart prior for InvQ
              delta = 2, # prior for Dirichlet distribution of w
              a=1,b=1, # parameters for gamma prior of tau, a = shape, b = rate
              aran=1,bran=1, # parameters for gamma prior of rantau, a = shape, b = rate
              V = diag(100, totnran), # Wishart prior scale matrix for InvQ
              ranmu0 = rep(0, totnran)) # prior mean for mu (mean of random effects)
names(param$ranmu0) <- rannames
param$gamma1 <- rep(-1, length(Ords)) # first fixed thresholds (one for each Ord outcome)
names(param$gamma1) <- Ords
param$gammaBin <- gamBin # one threshold for Binary outcomes (the same for all Bins)

param$fixmu <- param$fixD <- list()
for(y in c(Nums, Ords, Bins)){
  param$fixmu[[y]] <- rep(0, nfix[y]) # prior mean value of fixed effects
  param$fixD[[y]] <- rep(10, nfix[y]) # prior variances of fixed effects (later also multiplied by tau^{-1})
}
param$InvV <- backsolve(param$V, diag(1,totnran)) # Inverse of Wishart prior scale matrix for InvQ

# other used parameters (not for priors)
param$sd_b <- 0.01 # SD of random effects ONLY USED when sampling arbitrary initial values
param$sd_beta <- 0.01 # SD of fixed effects ONLY USED when sampling arbitrary initial values




###--------------------------------------------------------
### Initial values for sampling
###--------------------------------------------------------

# number of chains to be sampled 
Nchains <- 4

# Skip the following if you do not want to supply specific values from which to start
# the values below are based on the true parameter values from generating the dataset

### Inits based on true parameter values
# one list of initial values for each chain 
rinits <- list()
for(ch in 1:Nchains){
  rinits[[ch]] <- list()
  rinits[[ch]]$w <- rep(1/K, K)
  rinits[[ch]]$gamma <- list()
  rinits[[ch]]$gamma[[Ords[1]]] <- gam
  rinits[[ch]]$beta <- list()
  for(k in 1:K){
    rinits[[ch]]$beta[[k]] <- list()
    for(y in c(Nums, Ords, Bins)){
      rinits[[ch]]$beta[[k]][[y]] <- 
        Betamu[[K]][[k]][[r]][[d]][[gsub("f","",y)]][switch(r,
                                                            intercept = c(1,2,4),
                                                            slope = 1:3,
                                                            both = 1:2)]
    }
  }
  rinits[[ch]]$tau <- rep(tauNum, length(Nums))
  rinits[[ch]]$mu <- list()
  for(k in 1:K){
    if(r=="both"){
      rinits[[ch]]$mu[[k]] <- Eb[[K]][[k]][[r]][[d]][c(1,4,2,5,3,6)]
    }else{
      rinits[[ch]]$mu[[k]] <- Eb[[K]][[k]][[r]][[d]]
    }
  }
  rinits[[ch]]$rantau <- list()
  for(y in c(Nums, Ords, Bins)){
    rinits[[ch]]$rantau[[y]] <- rep(1, ifelse(r=="both",2,1))
  }
  if(spec["InvSigma"]){
    rinits[[ch]]$InvSigma <- list()
    for(k in 1:K){
      rinits[[ch]]$InvSigma[[k]] <- chol2inv(chol(Sigma[[K]][[k]][[r]]))
      if(r=="both"){
        rinits[[ch]]$InvSigma[[k]] <- chol2inv(chol(Sigma[[K]][[k]][[r]][c(1,4,2,5,3,6),c(1,4,2,5,3,6)]))
      }else{
        rinits[[ch]]$InvSigma[[k]] <- chol2inv(chol(Sigma[[K]][[k]][[r]]))
      }
    }
  }else{
    if(r=="both"){
      rinits[[ch]]$InvSigma <- chol2inv(chol(Sigma[[K]][[k]][[r]][c(1,4,2,5,3,6),c(1,4,2,5,3,6)]))
    }else{
      rinits[[ch]]$InvSigma <- chol2inv(chol(Sigma[[K]][[k]][[r]]))
    }
  }
  if(r=="both"){
    rinits[[ch]]$InvQ  <- Sigma[[K]][[k]][[r]][c(1,4,2,5,3,6),c(1,4,2,5,3,6)]*param$nu_0
  }else{
    rinits[[ch]]$InvQ  <- Sigma[[K]][[k]][[r]]*param$nu_0
  }
}

### Inits based on simulated dataset
for(ch in 1:Nchains){
  rinits[[ch]]$U <- data$U[seq(1, n*n_i, by = n_i)]-1
  rinits[[ch]]$pUik <- matrix(1/K, nrow = n, ncol = K)
  rinits[[ch]]$latent <- list()
  rinits[[ch]]$latent[["fY2"]] <- data$Y2
  rinits[[ch]]$latent[["Y3"]] <- data$Y3
  rinits[[ch]]$b <- data[seq(1,n*n_i,by=n_i), grep("b_Y", colnames(data))]
}



###---------------------------------------
### Gibbs sampling 
###---------------------------------------

B = 500   # length of burn-in period
M = 10000 # length of the final chain 
howsave = "matrix" # parameter describing the way how sampled parameters should be saved
# = "matrix" --> all parameters will be save in one big matrix under specific column names
# = "list" --> saved in structured lists

# Y - matrix of outcomes
# X - matrix of covariates
# K - number of clusters to be found 
# spec - what parameters should be cluster-specific
# calc - what additional parameters should be calculated
# whatsave - what parameters should be saved
# Nums - numeric outcomes 
# Ords - ordinal outcomes (0, 1, ..., L-1)
# Bins - binary outcomes (0 or 1)
# Id - name of column (both in Y and X) of subject id
# Formula - list of fixed/random formula - separately for each outcome (may be different)
# Nchains - number of chains to be sampled
# inits - list of Nchains initial values from which to start

## You do not need to specify spec, calc, whatsave or param.
## Function GibbsClassNumOrdBin has some default settings:
# cluster-specific: beta, mu, InvSigma   others common to all clusters
# calc all FALSE
# whatsave: F=latent,rantau,InvQ,b   others will be saved
# param: the same as used above

## Initial values supplied
# sampling Nchains chains - one after another (not in parallel)
# one can create a function that samples one chain with given inits
# to be called in parallel (see example below)
set.seed(123456)
system.time(
  mcmc <- GibbsClassNumOrdBin(Y = Y, X = X, K = K, 
                              spec = spec, calc = calc, whatsave = whatsave,
                              Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                              Formula = Formula,
                              M = M, B = B, Nchains = Nchains, 
                              inits = rinits, # initial values supplied
                              howsave = howsave)
)
# about 93 sec

## Continue in sampling, where the last chain stopped
# mcmc$last has the same structure as mcmc$inits
# it however contains values from the last step
set.seed(123456)
system.time(
  mcmc <- GibbsClassNumOrdBin(Y = Y, X = X, K = K, 
                              spec = spec, calc = calc, whatsave = whatsave,
                              Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                              Formula = Formula,
                              M = M, B = B, Nchains = Nchains, 
                              inits = mcmc$last, # mcmc$last as inits
                              howsave = howsave)
)
# about 93 sec

## OR DO IT IN PARALLEL
## function for parallel computation
SampleChainsParallel <- function(chain = 1, 
                                 n = n, Y = Y, XX = X, K = K,
                                 spec = spec, calc = calc, whatsave = whatsave,
                                 Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                                 Formula = Formula, inits = inits,
                                 M = M, B = B, howsave = howsave,
                                 RROOT, CROOT){
  chinits <- list()
  chinits[[1]] <- inits[[chain]]
  
  source(paste0(RROOT, "function_GibbsClassNumOrdBin.R"))
  source(paste0(RROOT, "FromCtoList.R"))
  source(paste0(RROOT, "FromListtoC.R"))
  source(paste0(RROOT, "FromCtoMatrix.R"))
  source(paste0(RROOT, "FromMatrixtoC.R"))
  
  library("mvtnorm")
  
  dyn.load(paste0(CROOT, "class_num_ord_bin.dll"))
  #dyn.load(paste0(CROOT, "class_num_ord_bin.so"))
  
  return(
    GibbsClassNumOrdBin(Y = Y, X = XX, K = K,
                        spec = spec, calc = calc, whatsave = whatsave,
                        Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                        Formula = Formula, inits = chinits,
                        M = M, B = B, Nchains = 1,
                        howsave = howsave)
  )
  
}

library(parallel)
myCluster <- makeCluster(Nchains) 
#myCluster <- makeCluster(Nchains, type = "FORK")      
clusterSetRNGStream(myCluster, 123456)

## parallel computation of chains
system.time(
  mcmcs <- parLapply(myCluster, 1:Nchains, SampleChainsParallel, 
                     n = n, Y = Y, XX = X, K = K,
                     spec = spec, calc = calc, whatsave = whatsave,
                     Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                     Formula = Formula, inits = mcmc$last,
                     M = M, B = B, howsave = howsave,
                     RROOT = RROOT, CROOT = CROOT)
)
# about 28 sec
stopCluster(myCluster)

### Getting the chains together
if(howsave == "matrix"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmcs[[chain]]$all[,"chain"] <- rep(chain, dim(mcmcs[[chain]]$all)[1])
      mcmc$all <- rbind(mcmc$all, mcmcs[[chain]]$all)
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}

if(howsave == "list"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmc[[chain]] <- mcmcs[[chain]][[1]]
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}



## What if I would like to start the sampling, 
## however, this time I would like Sigma to be cluster-specific?

# first change the cluster-specification
spec["InvSigma"] <- T

# and adjust the initial values
Sginits <- rinits
# the inits for InvSigma have to be changed to list of K matrices
for(ch in 1:Nchains){
  # initiate with the same matrix in all clusters
  IS <- Sginits[[ch]]$InvSigma
  Sginits[[ch]]$InvSigma <- list()
  for(k in 1:K){
    Sginits[[ch]]$InvSigma[[k]] <- IS
  }
}

set.seed(123456)
system.time(
  mcmcSg <- GibbsClassNumOrdBin(Y = Y, X = X, K = K, 
                                spec = spec, calc = calc, whatsave = whatsave,
                                Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                                Formula = Formula,
                                M = M, B = B, Nchains = Nchains, 
                                inits = Sginits, # mcmc$last as inits
                                howsave = howsave)
)

## Our generated data have actually beta common to all clusters
## but not in the model specification
# change the specificity of beta
spec["beta"] <- F
# also change back the cluster-specificity of InvSigma (due to example above)
spec["InvSigma"] <- F
# this time save allocations U
whatsave["U"] <- T

# and adjust the initial values
binits <- rinits
binits[[1]]$beta
# the inits for beta have to be changed from list of lists to just a list
for(ch in 1:Nchains){
  # initiate with the same matrix in all clusters
  binits[[ch]]$beta <- list()
  for(y in c(Nums, Ords, Bins)){
    k = 1 # with r=d=intercept should be the same for any k
    binits[[ch]]$beta[[y]] <- 
      Betamu[[K]][[k]][[r]][[d]][[gsub("f","",y)]][c(1,2,4)]
  }
}

set.seed(123456)
system.time(
  mcmcb <- GibbsClassNumOrdBin(Y = Y, X = X, K = K, 
                                spec = spec, calc = calc, whatsave = whatsave,
                                Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                                Formula = Formula,
                                M = M, B = B, Nchains = Nchains, 
                                inits = binits, # mcmc$last as inits
                                howsave = howsave)
)

## Initial values not supplied
# InitType is a vector of length Nchains declaring how it should be initialized
# (each chain can have different initialization rule)
# InitType == 1 <--> use linear regression to guess coefficients + estimate gammas
#                    (default)
# InitType == 2 <--> put zeros + estimate gammas
# InitType == 3 <--> generate something around zero + arithmetic sequence for gamma

# regardless of InitType the initial cluster allocations are sampled from Unif(0, 1, ..., K-1)
spec["beta"] <- T
set.seed(123456)
system.time(
  mcmc <- GibbsClassNumOrdBin(Y = Y, X = X, K = K, 
                              spec = spec, calc = calc, whatsave = whatsave,
                              Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                              Formula = Formula,
                              M = M, B = B, Nchains = Nchains, 
                              InitType = c(1,2,3,1), # declaring type of inicialization
                              howsave = howsave)
)
# see, different initial values
mcmc$inits[[1]]$beta
mcmc$inits[[2]]$beta
mcmc$inits[[3]]$beta
mcmc$inits[[4]]$beta

# last values
mcmc$last[[1]]$beta
mcmc$last[[2]]$beta
mcmc$last[[3]]$beta
mcmc$last[[4]]$beta

## What else does the output mcmc contain?
mcmc$howsave
mcmc$Nchains
mcmc$all[1:10, 1:10] # matrix of all sampled values (since howsave = "matrix")
dim(mcmc$all) # contains all four chains stacked below 42000 = (B+M)*Nchains
# (k) in the name means that this parameter is specific to cluster (k)
# [j] means that it is j-th element of that parameter
# both (k)[j] can be combined!

# description of model parameters
# cluster-specificity, outcome-specificity, dimensions
mcmc$settings
# used for creating column names

mcmc$B # burn-in period... states 1,...,B will not be used in plots 
mcmc$M # length of the MCMC sample
mcmc$BM # B+M
mcmc$Nums; mcmc$Ords; mcmc$Bins
mcmc$spec; mcmc$calc; mcmc$whatsave
mcmc$K
mcmc$param
mcmc$ncat # number of categories in ordinal outcomes
mcmc$Formula
# fixed effects (# by outcomes, total, names, names in list)
mcmc$nfix; mcmc$totnfix; mcmc$fixnames; mcmc$lfixnames
# random effects (# by outcomes, total, names, names in list)
mcmc$nran; mcmc$totnran; mcmc$rannames; mcmc$lrannames



###------------------------------------------------
### Plotting the output
###------------------------------------------------

## File "plotting_functions.R" contains several functions for plotting the output
## Implementation works with both types of saving the outputs
# plot.traceplots - samled values vs. iteration number
# plot.kerneldensity - kernel density estimate for given parameter
# plot.ECDF - empirical distribution function
# plot.ACF - auto-correlation function
# plot.all - 2×2 plot containing the four plots above
# plot.classes - compares values of cluster-specific parameter among clusters
#              - by either kernel density or ECDF

## Thinnig of sampled chains - to reduce the number of plotted values
figthin <- max(round(M/1000), 1)
ChainsToPlot <- 1:Nchains

## probabilites w
# w is not cluster specific, but rather a vector of length K
# hence: what = "w", dispec = "k", rest missing
pdf(paste0(FIG, "mcmc_w_K_", K, ".pdf"), 
    width = 7, height = 6)
{
  par(mfrow = c(K,4))
  for(k in 1:K){
    plot.traceplots(mcmc = mcmc, thin = figthin, 
                    what = "w", dimspec = k, labcex = 0.6, 
                    whatLAB = paste0("w[",k,"]"),
                    ChainsToPlot = ChainsToPlot)
    plot.kerneldensity(mcmc = mcmc, thin = figthin, 
                       what = "w", dimspec = k, labcex = 0.6, 
                       whatLAB = paste0("w[",k,"]"),
                       ChainsToPlot = ChainsToPlot)
    plot.ECDF(mcmc = mcmc, thin = figthin, 
              what = "w", dimspec = k, labcex = 0.6, 
              whatLAB = paste0("w[",k,"]"),
              ChainsToPlot = ChainsToPlot)
    plot.ACF(mcmc = mcmc, thin = figthin, 
             what = "w", dimspec = k, labcex = 0.6, 
             whatLAB = paste0("w[",k,"]"), MaxLag = 30,
             ChainsToPlot = ChainsToPlot)
  }
}
dev.off()

## sdSigma - standard deviations of Sigma matrix (square root of diag)
# not cluster-specific (mcmc$spec["InvSigma"] = F)
# vector of length totnran
# hence: what = "sdSigma", dimspec = j, where j in 1:totnran
totnran # dimension of matrix Sigma
# the following order of elements
rannames
Nums
Ords
Bins

pdf(paste0(FIG, "mcmc_sdSigma_K_", K, ".pdf"), 
    width = 7, height = ifelse(totnran == 3, 6, 9))
{
  par(mfrow = c(totnran,4))
  for(j in 1:totnran){
    plot.traceplots(mcmc = mcmc, thin = figthin, 
                    what = "sdSigma", dimspec = j, 
                    labcex = 0.6,
                    ChainsToPlot = ChainsToPlot)
    plot.kerneldensity(mcmc = mcmc, thin = figthin, 
                       what = "sdSigma", dimspec = j, 
                       labcex = 0.6,
                       ChainsToPlot = ChainsToPlot)
    plot.ECDF(mcmc = mcmc, thin = figthin, 
              what = "sdSigma", dimspec = j, 
              labcex = 0.6,
              ChainsToPlot = ChainsToPlot)
    plot.ACF(mcmc = mcmc, thin = figthin, 
             what = "sdSigma", dimspec = j, 
             labcex = 0.6, MaxLag = 30,
             ChainsToPlot = ChainsToPlot)
  }
}
dev.off()

## sdSigma from with cluster-specific Sigma matrix
# add: kspec = k, where k in 1:K
# and use mcmcSg which is a MCMC sample with Sigma cluster-specific
# (k) denotes cluster
# [j] denotes the element of sdSigma, see rannames
{
  par(mfrow = c(totnran,K))
  for(j in 1:totnran){
    for(k in 1:K){
      plot.traceplots(mcmc = mcmcSg, thin = figthin, 
                      what = "sdSigma", dimspec = j, kspec = k,
                      labcex = 0.6,
                      ChainsToPlot = ChainsToPlot)
    }
  }
}

## corSigma
# not cluster-specific
# symmetric non-diagonal matrix parameter (upper triangle)
mcmc$settings["corSigma",]
# hence: what = "corSigma", dimspec = c(i,j) with i <  j
pdf(paste0(FIG, "mcmc_corSigma_K_", K, ".pdf"), 
    width = 7, height = 6)
{
  par(mfrow = c(totnran-1, totnran-1))
  for(i in 1:(totnran-1)){
    for(j in 2:(totnran)){
      if(i < j){
        plot.traceplots(mcmc = mcmc, thin = figthin, 
                        what = "corSigma", dimspec = c(i,j), 
                        labcex = 0.6,
                        ChainsToPlot = ChainsToPlot)
      }else{
        plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "", 
             xaxt = "n", yaxt = "n", bty = "n")
      }
    }
  }
}
dev.off()

## Sigma
# not cluster-specific
# symmetric matrix with diagonal (upper triangle)
mcmc$settings["Sigma",]
# hence: what = "Sigma", dimspec = c(i,j) with i <=  j

pdf(paste0(FIG, "mcmc_Sigma_K_", K, ".pdf"), 
    width = 7, height = 6)
{
  par(mfrow = c(totnran,totnran))
  for (i in 1:totnran){
    for(j in (1:totnran)){
      if(i <= j){
        plot.traceplots(mcmc = mcmc, thin = figthin, 
                        what = "Sigma", dimspec = c(i,j), 
                        labcex = ifelse(totnran ==3, 0.6, 0.5),
                        ChainsToPlot = ChainsToPlot)
      }else{
        pomdata <- FindAllData(mcmc=mcmc, thin = figthin,
                               what = "Sigma", dimspec = c(j,i),
                               ChainsToPlot = ChainsToPlot)
        plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "", 
             xaxt = "n", yaxt = "n", bty = "n")
        text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
      }
    }
  }
}
dev.off()

## InvSigma
# not cluster-specific
# symmetric matrix with diagonal (upper triangle)
mcmc$settings["InvSigma",]
# hence: what = "InvSigma", dimspec = c(i,j) with i <=  j

pdf(paste0(FIG, "mcmc_InvSigma_K_", K, ".pdf"), 
    width = 7, height = 6)
{
  par(mfrow = c(totnran,totnran))
  for (i in 1:totnran){
    for(j in (1:totnran)){
      if(i <= j){
        plot.traceplots(mcmc = mcmc, thin = figthin, 
                        what = "InvSigma", dimspec = c(i,j), 
                        labcex = ifelse(totnran ==3, 0.6, 0.5),
                        ChainsToPlot = ChainsToPlot)
      }else{
        pomdata <- FindAllData(mcmc=mcmc, thin = figthin,
                               what = "InvSigma", dimspec = c(j,i),
                               ChainsToPlot = ChainsToPlot)
        plot(x=c(0,1), y=c(0,1), type = "n", xlab = "", ylab = "", 
             xaxt = "n", yaxt = "n", bty = "n")
        text(0.5,0.5, labels = format(mean(pomdata$y), digits = 2, nsmall = 2))
      }
    }
  }
}
dev.off()

## beta
# cluster-specific and outcome-specific vector parameter
# hence: what = "beta", kspec = k, yspec = y, dimspec = i
# (k) denotes cluster it belongs to
lfixnames
# [1] = effect of X1
# [2] = effect of X2
# [3] = effect of time
for(k in 1:K){
  # traceplots
  pdf(paste0(FIG, "traceplots_beta_K_", K, "_", k, ".pdf"), 
      width = 7, height = 6)
  {
    par(mfrow = c(sum(nY), nfix[1]))
    for (y in c(Nums, Ords, Bins)){
      for(i in 1:nfix[y]){
        plot.traceplots(mcmc = mcmc, thin = figthin, what = "beta", 
                        kspec = k, yspec = y, dimspec = i, labcex = 0.6,
                        ChainsToPlot = ChainsToPlot)
      }
    }
  }
  dev.off()
}
###!!! Note that each chain may have different interpretation of cluster lables!!!
### which happens when no inits are supplied and initial values are randomly sampled
### later we show how to deal with it

# beta - classes comparison - for each chain separately
# plot.classes by default plot both kernel and ECDF
# supress ECDF by doECDF = F
# supress parmfrow (default is c(1,2) for kernel and ECDF next to each other)

for(ch in 1:Nchains){
  pdf(paste0(FIG, "classes_beta_K_", K, "_ch_", ch, ".pdf"), 
      width = 7, height = 6)
  {
    par(mfrow = c(sum(nY), nfix[1]), mar = c(3.5,3.5,0.8,0.8))
    for(y in c(Nums, Ords, Bins)){
      for(i in 1:nfix[y]){
        plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, 
                     ChainsToPlot = ch,
                     what = "beta", yspec = y, dimspec = i, parmfrow = F)
      }
    }
  }
  dev.off()
}

## beta - version with not group-specific beta (use the other sample mcmcb)
# remove: kspec = k
{
  par(mfrow = c(sum(nY), nfix[1]))
  for (y in c(Nums, Ords, Bins)){
    for(i in 1:nfix[y]){
      plot.traceplots(mcmc = mcmcb, thin = figthin, what = "beta", 
                      yspec = y, dimspec = i, labcex = 0.6,
                      ChainsToPlot = ChainsToPlot)
    }
  }
}

{
  par(mfrow = c(sum(nY), nfix[1]))
  for (y in c(Nums, Ords, Bins)){
    for(i in 1:nfix[y]){
      plot.kerneldensity(mcmc = mcmcb, thin = figthin, what = "beta",
                         yspec = y, dimspec = i, labcex = 0.6,
                         ChainsToPlot = ChainsToPlot)
    }
  }
}


## mu
# cluster-specific vector parameter (but not outcome-specific like beta)
mcmc$settings["mu",]
# hence: what = "mu", kspec = k, dimspec = i and yspec missing
# (k) denotes cluster it belongs to
rannames
# random intercept terms for Nums, Ords, Bins - in that order
for(k in 1:K){
  # traceplots
  pdf(paste0(FIG, "traceplots_mu_K_", K, "_", k, ".pdf"), 
      width = 7, height = 6)
  {
    par(mfrow = c(totnran,4))
    for(j in 1:totnran){
      plot.traceplots(mcmc = mcmc, thin = figthin, 
                      what = "mu", kspec = k, dimspec = j, 
                      labcex = 0.6,
                      ChainsToPlot = ChainsToPlot)
      plot.kerneldensity(mcmc = mcmc, thin = figthin, 
                         what = "mu", kspec = k, dimspec = j, 
                         labcex = 0.6,
                         ChainsToPlot = ChainsToPlot)
      plot.ECDF(mcmc = mcmc, thin = figthin, 
                what = "mu", kspec = k, dimspec = j, 
                labcex = 0.6,
                ChainsToPlot = ChainsToPlot)
      plot.ACF(mcmc = mcmc, thin = figthin, 
               what = "mu", kspec = k, dimspec = j, 
               labcex = 0.6, MaxLag = 30,
               ChainsToPlot = ChainsToPlot)
    }
  }
  dev.off()
}
###!!! Note that each chain may have different interpretation of cluster labels!!!
### which happens when no inits are supplied and initial values are randomly sampled
### later we show how to deal with it

# mu - classes comparison - for each chain separately
# plot.classes by default plot both kernel and ECDF
# supress ECDF by doECDF = F
# supress parmfrow (default is c(1,2) for kernel and ECDF next to each other)

for(ch in 1:Nchains){
  pdf(paste0(FIG, "classes_mu_K_", K, "_ch_", ch, ".pdf"), 
      width = 7, height = 6)
  {
    par(mfrow = c(1, totnran), mar = c(3.5,3.5,0.8,0.8))
      for(i in 1:totnran){
        plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, 
                     ChainsToPlot = ch,
                     what = "mu", dimspec = i, parmfrow = F)
      }
  }
  dev.off()
}


## tau
# not cluster-specific, but Num-outcome-specific, only one-dimensional
# hence: what = "tau", yspec = Nums[1], else missing
pdf(paste0(FIG, "mcmc_K_", K, "_tau.pdf"), 
    width = 7, height = 6)
plot.all(mcmc = mcmc, thin = figthin,
         what = "tau", yspec = Nums[1], MaxLag = 30,
         ChainsToPlot = ChainsToPlot)
dev.off()

## gamma
# not cluster-specific, but Ord-outcome-specific, possibly vector
# hence: what = "gamma", yspec = Ords[1], dimspec = 1
y = Ords[1]
for(y in c(Ords)){
  pdf(paste0(FIG, "mcmc_gamma_K_", K, "_", y, ".pdf"), 
      width = 7, height = 6)
  plot.all(mcmc = mcmc, thin = figthin,
           what = "gamma", yspec = y, dimspec = 1, MaxLag = 30,
           ChainsToPlot = ChainsToPlot)
  dev.off()
}
# this parameter is the most problematic with respect to convergence
# hence, when continuing in sampling adjust the initial value appropriately


## If any of these traceplots shows that the convergence has not been met yet,
## then continue in sampling:
system.time(
  mcmc <- GibbsClassNumOrdBin(Y = Y, X = X, K = K, 
                              spec = spec, calc = calc, whatsave = whatsave,
                              Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                              Formula = Formula,
                              M = M, B = B, Nchains = Nchains, 
                              inits = mcmc$last, # mcmc$last as inits
                              howsave = howsave)
)

## OR IN PARALLEL
myCluster <- makeCluster(Nchains) 
#myCluster <- makeCluster(Nchains, type = "FORK")      
clusterSetRNGStream(myCluster, 123456)

system.time(
  mcmcs <- parLapply(myCluster, 1:Nchains, SampleChainsParallel, 
                     n = n, Y = Y, XX = X, K = K,
                     spec = spec, calc = calc, whatsave = whatsave,
                     Nums = Nums, Ords = Ords, Bins = Bins, Id = Id,
                     Formula = Formula, inits = mcmc$last,
                     M = M, B = B, howsave = howsave,
                     RROOT = RROOT, CROOT = CROOT)
)
# about 28 sec
stopCluster(myCluster)

### Getting the chains together
if(howsave == "matrix"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmcs[[chain]]$all[,"chain"] <- rep(chain, dim(mcmcs[[chain]]$all)[1])
      mcmc$all <- rbind(mcmc$all, mcmcs[[chain]]$all)
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}

if(howsave == "list"){
  mcmc <- mcmcs[[1]]
  if(Nchains > 1){
    for(chain in 2:Nchains){
      mcmc[[chain]] <- mcmcs[[chain]][[1]]
      mcmc$inits[[chain]] <- mcmcs[[chain]]$inits[[1]]
      mcmc$last[[chain]] <- mcmcs[[chain]]$last[[1]]
    }
  }
  mcmc$Nchains <- Nchains
}




###------------------------------------------------
### Relabelling cluster labels for whole chains
###------------------------------------------------

# When initializing with randomly generated partition
# the resulting cluster labels may lead to different permutations
# compare the changing colors interpretations in plotted classes for mu parameter
# among different chains
par(mfrow = c(2, 2), mar = c(3.5,3.5,0.8,0.8))
for(ch in 1:Nchains){
  for(i in 1:1){
    plot.classes(mcmc = mcmc, thin = figthin, doECDF = F, 
                 ChainsToPlot = ch,
                 what = "mu", dimspec = i, parmfrow = F)
  }
}

# Reorder so that 
# Red=1, Green=2, Blue=3, since even the real cluster labels are ordered in such a way
perm <- list()
perm[[1]] <- c(1, 3, 2)
perm[[2]] <- c(1, 2, 3)
perm[[3]] <- c(2, 3, 1)
perm[[4]] <- c(1, 2, 3)
# perm[[ch]][1] contains which old cluster should be the new first one (for chain ch)
# ... and so on

permmcmc <- PermuteClusterLabels(mcmc, perm)
# wait a while for relabelling...

par(mfrow = c(2, 2), mar = c(3.5,3.5,0.8,0.8))
for(ch in 1:Nchains){
  for(i in 1:1){
    plot.classes(mcmc = permmcmc, thin = figthin, doECDF = F, 
                 ChainsToPlot = ch,
                 what = "mu", dimspec = i, parmfrow = F)
  }
}

#!! All of this is possible only under NO label-switching problem!!
# We have the experience, that the chain "chooses" some interpretation
# of cluster labels and sticks with it. Hence, no label-switching occurs.
# However, it is not guaranteed that it will always happen.
# To check for label switching we use later method by Stevens.
# But first, clustering probabilities for all subjects and all iterations have to be evaluated...



###------------------------------------------------
### Coda compatibility 
###------------------------------------------------

#!! Perform only when NO label-switching occurs
#!! We also recommend to unify the cluster labels among different chains first (previous chunk)

### We have implemented a function that transfers our output
### into more familiar coda environment
codamcmc <- ToCodaMCMC(permmcmc, thin = 1)
# Especially useful for obtaining the summary statistics and quantiles
sumcodamcmc <- summary(codamcmc)
sumcodamcmc$statistics[grep("mu", rownames(sumcodamcmc$statistics)),]
sumcodamcmc$quantiles[grep("mu", rownames(sumcodamcmc$quantiles)),]

# maybe better separately for each chain
ch <- 1
sumch <- summary(codamcmc[[ch]])
sumch$statistics[grep("mu", rownames(sumch$statistics)),]
sumch$quantiles[grep("mu", rownames(sumch$quantiles)),]

# or plotting, however, commands below plot ALL parameters 
# and we do have a lot of them...
#traceplot(codamcmc)
#densplot(codamcmc)

###------------------------------------------------
### Cluster probabilities 
###------------------------------------------------

### An easy way for obtaining cluster probabilities
### is by the sampled U values
# remember, they have values in 0, ..., K-1
#Us <- mcmc$all[,c("chain", paste0("U[",1:n,"]"))] # original sample
Us <- permmcmc$all[,c("chain", paste0("U[",1:n,"]"))] # already permuted sample
clusteringU <- matrix(-1, nrow = as.numeric(n), ncol = Nchains)
certaintyU <- matrix(-1, nrow = as.numeric(n), ncol = Nchains)
for(i in 1:n){
  TAB <- table(Us[,i+1]+1, Us$chain)
  clusteringU[i,] <- as.numeric(rownames(TAB)[apply(TAB, 2, which.max)])
  certaintyU[i,] <- apply(TAB/mcmc$BM, 2, max)
}
# one column per each chain
# row i corresponds to subject i=1,...,n
# If ids are not directly 1,...,n then labelling subjects follows this permutation
UniqSubj <- unique(Y[,Id])
n <- length(UniqSubj)
NumUniqSubj <- c(1:n)
names(NumUniqSubj) <- UniqSubj

# clusteringU contains the cluster label, that has the highest frequency among BM states
head(clusteringU)

# certaintyU contains the relative frequency of the most frequent label
head(certaintyU)

for(ch in 1:Nchains){
  data[,paste0("clustering_",ch)] <- clusteringU[NumUniqSubj[data[,Id]],ch]
  data[,paste0("certainty_",ch)] <- certaintyU[NumUniqSubj[data[,Id]],ch]
}

# make the clustering 0, when certainty not large enough
for(ch in 1:Nchains){
  data[data[,paste0("certainty_",ch)] < 0.6, paste0("clustering_",ch)] <- 0
}

# compare it with the reality (use only the first observation per subject)
table(data[data$j==1,"clustering_1"], data[data$j==1,"U"])
table(data[data$j==1,"clustering_2"], data[data$j==1,"U"])
table(data[data$j==1,"clustering_3"], data[data$j==1,"U"])
table(data[data$j==1,"clustering_4"], data[data$j==1,"U"])
# meaning of the cluster labels may change from chain to chain and even differ from the reality
# this is caused by random initialization

### Approximation posterior distribution of classification probabilities
### by evaluation of p(U=k | theta) where latent variables are integrated out
### this involves integration of multivariate normal density over an interval
### which is solved by pmvnorm from library(mvtnorm)

# First we need to create a new dataset Y and a new dataset X
# containing data of subjects for which we want to evaluate the probabilities
# may be even data not used for creation of the model

# Since the calculation takes some time, we limit ourselves with
# only first nnew subjects and very large thinning
nnew <- 1000
Ynew <- Y[Y[,Id] <= nnew, ]     # outcomes 
Xnew <- X[X[,Id] <= nnew, ]     # covariates
Idnew <- Id   # label of column containing subject ids
probthin <- 200 # Only 10000/200 = 50 probabilities per subject evaluated

## Supply:
# mcmc - generated samples (where everything else needed is stored)
# Ynew, Xnew, Idnew
# thin
# start = mcmc$B+1, from which iteration number to start
# end = mcmc$BM, in which iteration number to end

## !! NOT RUN !! ##
# approximates chains one by one and not in parallel --> takes a lot of time
# system.time(
#   pmcmc <- pUnewk_pmvnorm(mcmc,
#                           Ynew = Ynew,
#                           Xnew = Xnew,
#                           Idnew = Idnew,
#                           thin = probthin)
# )

calculate_probs_for_chain <- function(ch, mcmc,
                                      Ynew, Xnew, Idnew,
                                      thin, start, end,
                                      RROOT, CROOT){
  source(paste0(RROOT, "FromCtoList.R"))
  source(paste0(RROOT, "FromListtoC.R"))
  source(paste0(RROOT, "FromCtoMatrix.R"))
  source(paste0(RROOT, "FromMatrixtoC.R"))
  source(paste0(RROOT, "pUnewk_C_pmvnorm.R"))
  dyn.load(paste0(CROOT, "pUnewk_pmvnorm.dll"))
  #dyn.load(paste0(CROOT, "pUnewk_pmvnorm.so"))
  
  library("mvtnorm")
  
  ## making mcmc as only 1 chain output
  chmcmc <- mcmc
  chmcmc$all <- mcmc$all[mcmc$all$chain == ch, ]
  chmcmc$all$chain <- 1
  chmcmc$Nchains <- 1
  
  RET <- pUnewk_pmvnorm(chmcmc,
                        Ynew = Ynew, Xnew = Xnew, Idnew = Idnew,
                        thin = thin, start = start, end = end)
  
  return(RET)
}

library(parallel)
myCluster <- makeCluster(mcmc$Nchains)
system.time(
  pmcmcs <- parLapply(myCluster, 1:mcmc$Nchains, calculate_probs_for_chain,
                      mcmc = permmcmc,
                      Ynew = Ynew, Xnew = Xnew, Idnew = Idnew,
                      thin = probthin, start = mcmc$B+1, end = mcmc$BM,
                      RROOT = RROOT, CROOT = CROOT)
)
# about 1200 sec for nnew=1000, probthin=200
stopCluster(myCluster)
colnames(pmcmcs[[1]])
# pUnewk - approximated probability
# pUerror - error of approximation from pmvnorm
# pUinform - inform indicator from pmvnorm (0 = converged, 1 = problems)

COL <- rainbow_hcl(K, c = 80, l = 50)
ch <- 1 # choose some chain

# plotting the posterior estimates
par(mfrow = c(2,5))
for(i in 1:10){
  plot(0,0,type = "n", xlim = c(0,1), ylim = c(0,10))
  for(k in 1:K){
    lines(density(pmcmcs[[ch]][,paste0("pUnewk[",i,",",k,"]")]),
          col = COL[k])
  }
}

### Classification based on HPD
# uses function HPDinterval from library(coda)
clusterHPD <- list()
for(ch in 1:mcmc$Nchains){
  chp <- pmcmcs[[ch]][, grep("pUnewk", colnames(pmcmcs[[ch]]))]
  clusterHPD[[ch]] <- numeric(nnew)
  
  for(i in 1:nnew){
    ichp <- chp[,paste0("pUnewk[",i,",",1:mcmc$K,"]")]
    
    HPDs <- HPDinterval(as.mcmc(ichp), prob = 0.6)
    means <- apply(ichp, 2, mean, na.rm = T)
    # cluster with maximal probability
    maxcl <- which.max(means)
    # however, it will be assigned there only if corresponding lower bound
    # exceeds any other upper bound
    cli <- ifelse(HPDs[maxcl, "lower"] > max(HPDs[setdiff(1:mcmc$K, maxcl),"upper"]),
                  maxcl, # lower bound of maxcl is higher than upper bound of others
                  0 # lower and upper bounds do cross --> not to be clustered
    )
    clusterHPD[[ch]][i] <- cli
    
  }
}

# 0 means not clustered
clusterHPD[[1]]
clusterHPD[[2]]
clusterHPD[[3]]
clusterHPD[[4]]

# Comparison with the reality 
for(ch in 1:Nchains){
  print(table(
    data[data$j==1 & data$subject<=nnew, "U"],
    clusterHPD[[ch]]))
}

# Comparison with the clustering based on sampled U values
for(ch in 1:Nchains){
  print(table(
    data[data$j==1 & data$subject<=nnew, paste0("clustering_",ch)],
    clusterHPD[[ch]]))
}



### How to compute classification probabilities dynamically?
# use function EachIdGradually
# takes each subject and creates subset of first j observations, j = 1, ..., n_i
# newid = (subject id)*10^digits_added + j
# --> each subject creates n_i "new" subjects
dyndata <- EachIdGradually(data, Id = Id, Order = "time", digits_added = 2)
Ydyn <- dyndata[, c("Y1", "fY2", "fY3", "subject", "newid")]
Ydyn$Y3 <- as.numeric(as.character(Ydyn$fY3))
Xdyn <- dyndata[, c("X1", "X2", "intercept", "time", "subject", "newid")]
Iddyn <- "newid"
dynids <- unique(Ydyn$newid)
ndyn <- length(dynids)
newidto1ton <- c(1:ndyn)
names(newidto1ton) <- dynids

### Classification probabilities calculation
## !! NOT RUN !!
# pmcmc <- pUnewk_pmvnorm(mcmc,
#                         Ynew = Ydyn,
#                         Xnew = Xdyn,
#                         Idnew = Iddyn,
#                         thin = probthin)

# and then use newidto1ton to link 1:ndyn used in names of pmcmc columns
# with corresponding newid
# and then with the original one
# ...



###------------------------------------------------
### Label-switching problem
###------------------------------------------------

### Method of Stevens (2000) for dealing with label-switching 
### works with classification probabilities.
### Hence we would need to calculate them with nnew = n and probthin = 1.
### We will showcase the use of this only to the thinned chain.

# set the precision and maximum number of iterations
eps <- 1e-3
max_iter <- 10

# First of all - matrix of all permutations on K elements
permutations <- function(n){
  if(n==1){
    return(matrix(1))
  }else{
    sp <- permutations(n-1)
    p <- nrow(sp)
    A <- matrix(nrow=n*p,ncol=n)
    for(i in 1:n){
      A[(i-1)*p+1:p,] <- cbind(i,sp+(sp>=i))
    }
    return(A)
  }
}
perms <- matrix(c(1:K)[permutations(K)],ncol=K)
perms <- perms - 1 # to be from 0 to K-1 and not 1,...,K

# select chain
chain <- 1
prob <- pmcmcs[[ch]][,grep("pUnewk", colnames(pmcmcs[[ch]]))]
summary(unlist(prob))
# if some NA, then replace them with 1/K
prob[is.na(prob)] <- 1/K
dim(prob)
prob[1:5, 1:9]

# for the following C function reorder columns
prob_reorder <- rep(seq(1,nnew*K,by=nnew), nnew) + rep(0:(nnew-1), each = K)
prob <- prob[,prob_reorder]
prob[1:5, 1:9]
#now probs for one subject are next to each other

results <-
  .C("FindPermutations",
     K = as.integer(K),           # IN  [1]     number of clusters
     M = as.integer(dim(prob)[1]),# IN  [1]     number of iterations
     n = as.integer(nnew),           # IN  [1]     number of subjects
     nperms = as.integer(dim(perms)[1]),    
     # IN  [1]     K! - number of permutations
     perms = as.integer(c(perms)),     
     # IN  [K!*K]  rows contain all permutations 
     probs = as.double(c(as.matrix(prob))),     
     # IN  [M*K*n] generated probabilities
     iteration = as.integer(0), 
     # OUT [1]     number of iterations needed to converge
     max_iter = as.integer(max_iter), 
     # IN  [1]     maximal number of iterations allowed
     nu = as.integer(c(matrix(rep(0:(K-1), each = dim(prob)[1]), ncol = K))),        
     # OUT [M*K]   m-th row contains permutation of m-th iteration
     Q = as.double(c(matrix(0, nrow = K, ncol = n))),      
     # OUT [K*n]   i-th column contains Q probabilities of i-th subject
     eps = as.double(eps)    
     # IN  [1]     tolerance for convergence
  )

nu <- matrix(results$nu, nrow = dim(prob)[1], ncol = K)
dim(nu)
head(nu)
# each row contains a suitable permutation of cluster labels for corresponding iteration
# In our experience, no label-switching problem is met
# cluster labels are static for the chain during sampling
summary(nu)
# if values in the same column are the same --> no correction for label switching needed

# Otherwise you should permute all cluster-specific parameters 
# at each iteration by the corresponding row of nu
# not yet implemented .... since not needed, yet....



