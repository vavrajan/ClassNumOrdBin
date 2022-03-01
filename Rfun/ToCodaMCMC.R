### Function for transferring GibbsClassNumOrdBin output
### into output of library("coda") style

ToCodaMCMC <- function(mcmc, thin = 1){
  if(mcmc$howsave == "list"){
    # Transferring to shape of input as C
    params <- c("w", "U", "pUik", "pUik_nonb", "latent", "tau",  
                "mu", "InvSigma", "detInvSigma", "detSigma", "Sigma", 
                "sdSigma", "corSigma", "InvQ", "detInvQ", "detQ", "Q", "b")
    otherparam <- c("gamma", "min", "max", "rantau", "beta")
    
    cmcmc <- list()
    for(ch in 1:Nchains){
      cmcmc[[ch]] <- list()
      for(p in c(params,otherparam)){
        if(mcmc$settings[p, "save"]){
          if(is.element(p, otherparam)){
            cmcmc[[ch]][[p]] <- FromListtoC_settings(mcmc[[ch]], p=p, 
                                                     settings = mcmc$settings,
                                                     yspecd1 = mcmc$yspecd1[[p]])
          }else{
            cmcmc[[ch]][[p]] <- FromListtoC_settings(mcmc[[ch]], p=p, 
                                                     settings = mcmc$settings)
          }
        }
      } # end of for p in params and otherparams
    } # end of for ch in 1:Nchains
    
    # Transferring to shape of output as howsave=="matrix"
    newmcmc <- list()
    newmcmc$all <- data.frame()
    for(ch in 1:Nchains){
      AllData <- matrix(ch, nrow = B+M, ncol = 1)
      colnames(AllData) <- "chain"
      
      for(p in c(params,otherparam)){
        if(mcmc$settings[p, "save"]){
          if(is.element(p, otherparam)){
            AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[ch]][[p]],
                                                             p = p,
                                                             settings = mcmc$settings,
                                                             yspecd1 = mcmc$yspecd1[[p]]))
          }else{
            AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[ch]][[p]],
                                                             p = p,
                                                             settings = mcmc$settings))
          }
        }
      }
      
      newmcmc$all <- rbind(newmcmc$all, AllData)
    }
    newmcmc$B <- mcmc$B
    newmcmc$M <- mcmc$M
    newmcmc$BM <- mcmc$BM
    newmcmc$Nchains <- mcmc$Nchains
    mcmc <- newmcmc
  }
  
  # Now mcmc is of the shape of howsave == "matrix"
  codamcmc <- list()
  for(ch in 1:mcmc$Nchains){
    codamcmc[[ch]] <- mcmc$all[mcmc$all$chain == ch, 2:dim(mcmc$all)[2]]
    codamcmc[[ch]] <- mcmc(codamcmc[[ch]][seq(mcmc$B+1, mcmc$BM, by = thin),], thin = thin, start = mcmc$B+1, end = mcmc$BM)
  }
  codamcmc <- as.mcmc.list(codamcmc)
  return(codamcmc)
}
