### Function that takes new subjects, their outcomes and regressors
# also takes previously generated parameters from posterior distribution
# and calculates probabilities of classification

#dyn.load(paste(gsub("EU-SILC", "Classification", WD), "pUnewk_pmvnorm.dll", sep = "/"))


pUnewk_pmvnorm <- function(mcmc,      # MCMC generated parameters
                           #   output of GibbsClassNumOrdBin
                           Ynew,      # outcomes of new subjects
                           Xnew,      # regressors of new subjects
                           Idnew,     # Identification numbers
                           # common to rows of Ynew and Xnew 
                           # when they should be taken together
                           # Some subjects may be calculated gradually
                           # in such a case the Id is constructed as follows
                           # Id * 100 + number of rows
                           start = mcmc$B+1,
                           end = mcmc$BM,
                           thin = 1
                           )
{
  ### Every auxiliary values are stored within mcmc output
  # Such as whatsave, calc, spec, ...
  
  # Iterations to be taken
  iterations <- seq(start, end, by = thin)
  niter <- length(iterations)
  
  # For probability calculation we need the following parameters:
  #   w, tau, beta, mu, InvSigma
  # If those were not saved during the GibbsClassNumOrdBin function,
  # then we cannot use this function to calculate it.
  if(!(mcmc$whatsave["w"] & mcmc$whatsave["gamma"] & mcmc$whatsave["tau"] 
       & mcmc$whatsave["beta"] & mcmc$whatsave["mu"] & mcmc$whatsave["InvSigma"])){
    stop("Some of the parameters w, gamma, tau, beta, mu, InvSigma was not saved.
Change the whatsave settings and use GibbsClassNumBinOrd again.")
  }
  
  # Calculation will be done by C function, to speed up calculation
  # parameters need to be given to C function in a similar way to GibbsClassNumOrdBin
  
  nnew <- length(unique(Ynew[,Idnew]))
  Nnew <- dim(Ynew)[1]
  UniqSubjnew <- unique(Ynew[,Idnew])
  nnew <- length(UniqSubjnew)
  NumUniqSubjnew <- c(1:nnew)
  names(NumUniqSubjnew) <- UniqSubjnew
  nsubjnew <- table(Ynew[,Idnew])
  
  ### in case howsave == matrix, create place where to store results
  if(mcmc$howsave == "matrix"){
    output <- data.frame()
  }
  
  ### in case howsave = list, create a list of places to store results
  if(mcmc$howsave == "list"){
    output <- list()
  }
  
  ### Settings matrix
  params <- c("pUnewk", "pUerror", "pUinform", "w", "tau", "mu", "InvSigma")
  otherparam <- c("gamma", "beta")
  
  settings <- data.frame(save = sapply(c(params,otherparam), 
                                       function(p){
                                         if(is.element(p, c(names(mcmc$whatsave), names(mcmc$calc)))){
                                           return(c(na.omit(c(mcmc$whatsave[p], mcmc$calc[p]))))
                                         }else{
                                           return(T)
                                         }}),
                         isspec = sapply(c(params,otherparam), 
                                         function(p){
                                           if(is.element(p, c(names(mcmc$whatsave), names(mcmc$calc)))){
                                             return(ifelse(is.na(mcmc$spec[p]), F, mcmc$spec[p]))
                                           }else{
                                             return(F)
                                           }}),
                         K = mcmc$K,
                         isy = F,
                         ynums = 0,
                         yords = 0,
                         ybins = 0,
                         BM = mcmc$B+mcmc$M,
                         d1 = 0,
                         d2 = 0, 
                         BYROW = T,
                         sym = 0,
                         diag = 0,
                         diagval = 1,
                         D = 0)
  
  rownames(settings) <- c(params, otherparam)
  #settings
  
  # Individual changes
  settings["pUnewk",  c("BM", "d1", "d2", "BYROW", "D")] = c(niter, nnew, mcmc$K, F, 2)
  settings["pUerror", c("BM", "d1", "d2", "BYROW", "D")] = c(niter, nnew, mcmc$K, F, 2)
  settings["pUinform", c("BM", "d1", "d2", "BYROW", "D")] = c(niter, nnew, mcmc$K, F, 2)
  settings["w", c("d1", "D")] = c(mcmc$K, 1)
  settings["gamma", c("isy", "yords", "D")] = c(T,1,1)
  settings["tau", c("isy", "ynums")] = c(T, 1)
  settings["beta", c("isy", "ynums", "yords", "ybins", "D")] = c(T,1,1,1,1)
  settings["mu", c("d1", "D")] = c(mcmc$totnran, 1)
  settings["InvSigma", c("d1","d2", "sym", "diag", "D")] = c(mcmc$totnran,mcmc$totnran,1,1,2)
  
  ### Preparation of parameters for calculation of pUnewk in C
  ## Data
  # First column is going to be the id variable (0-th column)
  cId <- NumUniqSubjnew[as.character(Ynew[,Idnew])] - 1
  # -1 is there for C which works better with number beggining with 0
  # other columns  (beggining with 1st column)
  cY <- cX <- numeric()
  for(y in c(mcmc$Nums, mcmc$Ords, mcmc$Bins)){
    if(!is.element(y, colnames(Ynew))){
      stop(paste("You did not deliver outcome:", y))
    }
    cY <- c(cY, as.numeric(as.character(Ynew[,y])))
  }
  for(x in setdiff(colnames(Xnew), Idnew)){
    cX <- c(cX, Xnew[,x])
  } # no need for Id
  # Formula
  cFormulaF <- cFormulaR <- numeric()
  for(y in c(mcmc$Nums, mcmc$Ords, mcmc$Bins)){
    # regressors for fixed part
    for(i in 1:mcmc$nfix[y]){
      if(!is.element(mcmc$Formula[[y]]$fixed[i], colnames(Xnew))){
        stop(paste("You did not deliver regressor:", mcmc$Formula[[y]]$fixed[i]))
      }
    }
    cFormulaF <- c(cFormulaF, which(is.element(setdiff(colnames(Xnew), Idnew), mcmc$Formula[[y]]$fixed)))
    
    # regressors for random part
    for(i in 1:mcmc$nran[y]){
      if(!is.element(mcmc$Formula[[y]]$random[i], colnames(Xnew))){
        stop(paste("You did not deliver regressor:", mcmc$Formula[[y]]$random[i]))
      }
    }
    cFormulaR <- c(cFormulaR, which(is.element(setdiff(colnames(Xnew), Idnew), mcmc$Formula[[y]]$random)))
  }
  cFormulaF <- cFormulaF - 1
  cFormulaR <- cFormulaR - 1
  cnfix <- as.numeric(mcmc$nfix) # number of FIXED  regressors with variables y
  cnran <- as.numeric(mcmc$nran) # number of RANDOM regressors with variables y
  # dims - in the following order:
  # w, pUik, U can be determined out of K and n
  cdims <- c(mcmc$K, # w
             nnew*mcmc$K, # pUnewk
             sum(mcmc$ncat-2), # gamma
             mcmc$totnfix, # beta
             mcmc$nY["Nums"], # tau
             mcmc$totnran, # mu
             mcmc$totnran*(mcmc$totnran+1)/2 # InvSigma
  )
  cdimswithK <- c(mcmc$K, # w
                  nnew*mcmc$K, # pUnewk
                  sum(mcmc$ncat-2)*ifelse(mcmc$spec["gamma"], mcmc$K, 1), # gamma
                  mcmc$totnfix*ifelse(mcmc$spec["beta"], mcmc$K, 1), # beta
                  mcmc$nY["Nums"]*ifelse(mcmc$spec["tau"], mcmc$K, 1), # tau
                  mcmc$totnran*ifelse(mcmc$spec["mu"], mcmc$K, 1), # mu
                  mcmc$totnran*(mcmc$totnran+1)/2*ifelse(mcmc$spec["InvSigma"], mcmc$K, 1) # InvSigma
  )
  names(cdims) <- names(cdimswithK) <- c("w","pUnewk","gamma","beta","tau","mu","InvSigma")
  for(p in names(cdims)){
    settings[p, "dims"] = cdims[p]
    settings[p, "dimswithK"] = cdimswithK[p]
  }
  settings["pUerror", "dims"] <- settings["pUerror", "dimswithK"] <- 
    settings["pUinform", "dims"] <- settings["pUinform", "dimswithK"] <- nnew*mcmc$K
  # those dimensions are calculable in C function as well, when
  # K, n, ncat, nY, spec, nran, nfix     are available
  
  for(ch in 1:mcmc$Nchains){
    # arrays/fields where those variables are and will be stored
    cpUnewk     <- double(niter*cdimswithK[2])
    cpUerror    <- double(niter*cdimswithK[2])
    cpUinform   <- double(niter*cdimswithK[2])
    if(mcmc$howsave == "list"){
      cw        <- FromListtoC_settings(mcmc[[ch]], iterations, "w", settings)
      cgamma    <- FromListtoC_settings(mcmc[[ch]], iterations, "gamma", settings, yspecd1 = mcmc$ncat - 2)
      cbeta     <- FromListtoC_settings(mcmc[[ch]], iterations, "beta", settings, yspecd1 = mcmc$nfix)
      ctau      <- FromListtoC_settings(mcmc[[ch]], iterations, "tau", settings)
      cmu       <- FromListtoC_settings(mcmc[[ch]], iterations, "mu", settings)
      cInvSigma <- FromListtoC_settings(mcmc[[ch]], iterations, "InvSigma", settings)
    }
    if(mcmc$howsave == "matrix"){
      cw        <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "w", settings)
      cgamma    <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "gamma", settings, yspecd1 = mcmc$ncat - 2)
      cbeta     <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "beta", settings, yspecd1 = mcmc$nfix)
      ctau      <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "tau", settings)
      cmu       <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "mu", settings)
      cInvSigma <- FromMatrixtoC_settings(mcmc$all, chain = ch, iterations, "InvSigma", settings)
    }
    
    #dyn.load(paste(gsub("EU-SILC", "Classification", WD), "pUnewk_calculation.dll", sep = "/"))
    #dyn.unload(paste(gsub("EU-SILC", "Classification", WD), "pUnewk_calculation.dll", sep = "/"))
    
    #system.time(
    cmcmc <-
      .C("pUnewk_pmvnorm",
         Id        = as.integer(cId),
         Y         = as.double(cY),
         X         = as.double(cX),
         spec      = as.integer(mcmc$spec),
         # parameters describing dimensions
         chain     = as.integer(ch), # number of the chain
         K         = as.integer(mcmc$K), # number of classes
         BM        = as.integer(niter), # total number of generated states
         FormulaF  = as.integer(cFormulaF), # numbers of columns of X that should be used for FIXED  effects of modelled responses
         FormulaR  = as.integer(cFormulaR), # numbers of columns of X that should be used for RANDOM effects of modelled responses
         Nnew      = as.integer(Nnew), # total number of observations
         dims      = as.integer(cdims), # the length of subarray that corresponds to one state (disected by various parameters)
         dimswithK = as.integer(cdimswithK), # the length of subarray that corresponds to one state 
         # (disected by various parameters, also multiplication by K incorporated when such parameters is class-specific)
         ncat      = as.integer(mcmc$ncat), # the counts of categories of ordinal variables
         nY        = as.integer(mcmc$nY), # 3 numbers: counts of Nums, Ords and Bins variables
         nnew      = as.integer(nnew), # total number of subjects (different ids in the dataset)
         nfix      = as.integer(cnfix), # 
         nran      = as.integer(cnran),
         #predictor = double(N*sum(nY)),
         # the function should count totnran, totnfix and cumsum versions of nfix and nran
         # arrays to store generated states
         w           = as.double(cw),
         pUnewk      = as.double(cpUnewk),
         pUerror     = as.double(cpUerror),
         pUinform    = as.double(cpUinform),
         gamma       = as.double(cgamma),
         gamma1      = as.double(mcmc$param$gamma1),
         gammaBin    = as.double(mcmc$param$gammaBin),
         beta        = as.double(cbeta),
         tau         = as.double(ctau),
         mu          = as.double(cmu),
         InvSigma    = as.double(cInvSigma)
    )
    
    if(mcmc$howsave == "list"){
      # results will be returned in structured list - the same way as original R function
      chain <- list()
      for(p in c("pUnewk", "pUerror", "pUinform")){
        if(settings[p, "save"]){
          chain[[p]] <- FromCtoList_settings(cmcmc[[p]], p, settings)
        }
      }
      
      output[[ch]] <- chain
    }
    
    if(mcmc$howsave == "matrix"){
      # results will be returned in structured list - the same way as original R function
      AllData <- matrix(ch, nrow = niter, ncol = 1)
      colnames(AllData) <- "chain"
      
      for(p in c("pUnewk", "pUerror", "pUinform")){
        if(settings[p, "save"]){
          AllData <- cbind(AllData, FromCtoMatrix_settings(cmcmc[[p]], p, settings,
                                                           Nums = Nums, Ords = Ords, Bins = Bins))
        }
      }
      
      output <- rbind(output, AllData)
    }
    
  }  # end of for ch in 1:Nchains
  
  return(output)
  
}