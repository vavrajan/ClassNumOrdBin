#dyn.load(paste(gsub("EU-SILC", "Classification", WD), "class_num_ord_bin.dll", sep = "/"))
#source(paste(gsub("EU-SILC", "Classification", WD), "FromCtoList.R", sep = "/"))

### R Function of Gibbs sampling
GibbsClassNumOrdBin <- function(Y, X, K=2, spec, calc, whatsave,
                                Nums, Ords, Bins, Id, Formula,
                                M = 1000, B = 1000, Nchains = 1,
                                param, inits, InitType, howsave){
  mcmc <- list()
  N <- dim(Y)[1]
  UniqSubj <- unique(Y[,Id])
  n <- length(UniqSubj)
  NumUniqSubj <- c(1:n)
  names(NumUniqSubj) <- UniqSubj
  nsubj <- table(Y[,Id])
  
  ncat <- numeric(length(Ords))
  names(ncat) <- Ords
  lev <- list()
  for(o in Ords){
    ncat[o] <- nlevels(as.factor(Y[,o]))
    lev[[o]] <- levels(as.factor(Y[,o]))
  }
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
  totnran <- sum(nran)
  totnfix <- sum(nfix)
  
  if(!all(Y[,Id] == X[,Id])){stop("Rows of Y and X do not match.")}
  
  if (missing(spec)){
    spec <- c(F, F, T, F, T, F, T)
    names(spec) <- c("gamma", "tau", "beta", "rantau", "mu", "InvQ", "InvSigma")
  }
  if(spec["tau"]){spec["beta"] <- T}
  if(spec["rantau"]){spec["mu"] <- T}
  if(spec["InvQ"]){spec["InvSigma"] <- T}
  
  if(missing(calc)){
    calc <- c(F,
              F,F,F,
              F,F,F,F,F)
    names(calc) <- c("pUik_nonb", 
                     "Q", "detQ", "detInvQ",
                     "Sigma", "detSigma", "detInvSigma", "sdSigma", "corSigma")
  }
  
  if(missing(whatsave)){
    whatsave <- c(T,T,T,
                  T,T,T,F,
                  T,T,
                  T,F,
                  T,F,F)
    names(whatsave) <- c("w", "U", "pUik",
                         "gamma", "min", "max", "latent",
                         "beta", "tau", 
                         "mu", "rantau",
                         "InvSigma", "InvQ", "b")
  }
  
  if(!(whatsave["w"] & whatsave["gamma"] & whatsave["tau"] 
       & whatsave["beta"] & whatsave["mu"] & whatsave["InvSigma"])){
    warning("Some of the parameters w, gamma, tau, beta, mu, InvSigma will not be saved.
Therefore, you will not be able to use pUnewk_calculation for calculating probabilities of new subjects.
Change the whatsave settings and use GibbsClassNumBinOrd again if you want to calculate these probabilities afterwards.")
  }
  
  if(missing(howsave)){
    howsave = "list"
  }
  
  if (missing(param)){
    param <- list(nu_0 = totnran + 1,
                  nu_1 = totnran + 1,
                  delta = 2,
                  a=1,b=1,
                  aran=1,bran=1,
                  V = diag(100, totnran),
                  ranmu0 = rep(0, totnran))
    names(param$ranmu0) <- rannames
    param$gamma1 <- rep(-1, length(Ords))
    names(param$gamma1) <- Ords
    param$gammaBin <- 0
    param$fixmu <- param$fixD <- list()
    for(y in c(Nums, Ords, Bins)){
      param$fixmu[[y]] <- rep(0, nfix[y])
      param$fixD[[y]] <- rep(10, nfix[y])
    }
    param$InvV <- backsolve(param$V, diag(1,totnran))
    param$sd_b <- 0.01
    param$sd_beta <- 0.01
  } # end of if missing(param)
  
  ## Useful subdata and precalculated variables
  Xfix <- XtX <- Z <- tZ <- ZtZ <- XtXplusD <- cholXtXplusD <- list()
  for(y in c(Nums,Ords,Bins)){
    Xfix[[y]] <- as.matrix(X[, Formula[[y]]$fixed])
    XtX[[y]] <- t(Xfix[[y]]) %*% Xfix[[y]]
    Z[[y]] <- as.matrix(X[, Formula[[y]]$random])
    ZtZ[[y]] <- t(Z[[y]]) %*% Z[[y]]
    if(nfix[y]>1){
      XtXplusD[[y]] <- XtX[[y]]+diag(1/param$fixD[[y]])
      cholXtXplusD[[y]] <- chol(XtXplusD[[y]])
    }else{ 
      XtXplusD[[y]] <- XtX[[y]]+1/param$fixD[[y]]
      cholXtXplusD[[y]] <- sqrt(XtXplusD[[y]])
      colnames(cholXtXplusD[[y]]) <- rownames(cholXtXplusD[[y]]) <- Formula[[y]]$fixed
    }
  }
  for(j in 1:n){
    Z[[j]] <- tZ[[j]] <- ZtZ[[j]] <- list()
    for(y in c(Nums, Ords, Bins)){
      Z[[j]][[y]] <- as.matrix(X[Y[,Id] == UniqSubj[j],Formula[[y]]$random])
      colnames(Z[[j]][[y]]) <- Formula[[y]]$random
      tZ[[j]][[y]] <- t(Z[[j]][[y]])
      ZtZ[[j]][[y]] <- tZ[[j]][[y]] %*% Z[[j]][[y]]
    }
  }
  Ordasnum <- list()
  for(o in Ords){
    #Ordasnum[[o]] <- as.numeric(as.character(Y[,o]))
    Ordasnum[[o]] <- as.numeric(factor(Y[,o]))
    #if (all(is.na(Ordasnum[[o]]))){Ordasnum[[o]] <- as.numeric(Y[,o])}
  }
  ##
  
  initsgiven = 1
  ### Inicialization by rough estimates (chain 1)
  ### Inicialization by zeros (chain 2)
  ### Inicialization by almost zeros (chain 3)
  if (missing(inits)){
    cat("Inits are missing, calculating...\n")
    initsgiven = 0
    if (missing(InitType)){InitType <- rep(1, Nchains)}
    ## InitType == 1 <--> use linear regression to guess coefficients + estimate gammas
    ## InitType == 2 <--> put zeros + estimate gammas
    ## InitType == 3 <--> generate something around zero + arithmetic sequence for gamma
    inits <- list()
    for(ch in 1:Nchains){
      inits[[ch]] <- list()
      # parameters connected to classification
      inits[[ch]]$w <- rep(1/K, K)
      inits[[ch]]$U <- sapply(1:n, function(x){
        newU <- sump <- 0
        u1 <- runif(1)
        while(sump < u1){
          newU <- newU+1
          sump <- sump + inits[[ch]]$w[newU]
        }
        return(newU-1) # needs values in 0, ..., K-1
      })
      Nk <- selection <- list()
      for(k in 1:K){
        Nk[[k]] <- UniqSubj[inits[[ch]]$U+1 == k]
        selection[[k]] <- is.element(Y[,Id], Nk[[k]])
      }
      inits[[ch]]$pUik <- matrix(rep(inits[[ch]]$w, n), ncol = K, nrow = n, byrow = T)
      # parameters connected to latent variable modelling
      inits[[ch]]$gamma <- inits[[ch]]$latent <- list()
      if(spec["gamma"]){for(k in 1:K){inits[[ch]]$gamma[[k]] <- list()}}
      for(o in Ords){
        if (InitType[ch] == 3){
          # we do not have any latent data to actually be able to estimate roughly gamma
          # we will use an arithmetic sequence with step of size 1 (can be changed)
          if(spec["gamma"]){for(k in 1:K){inits[[ch]]$gamma[[k]][[o]] <- seq(param$gamma1[o]+1, param$gamma1[o]+ncat[o]-2, by = 1)}
          }else{inits[[ch]]$gamma[[o]] <- seq(param$gamma1[o]+1, param$gamma1[o]+ncat[o]-2, by = 1)}
          # we do not need initial values for latent variables of ordinal variables
          inits[[ch]]$latent[[o]] <- Ordasnum[[o]]+param$gamma1[o]-1.5
        }else{ # we will estimate gammas by the empirical probabilities
          if(spec["gamma"]){
            inits[[ch]]$latent[[o]] <- numeric(N)
            for(k in 1:K){
              probs <- cumsum(table(Y[selection[[k]],o]))
              probs <- probs/sum(probs)
              tildemu <- param$gamma1[o]-qnorm(probs[1])
              inits[[ch]]$gamma[[k]][[o]] <- tildemu + qnorm(probs[2:(ncat[o]-1)])
              if(ncat[o]>3){
                Vals <- c(param$gamma1[o] - 0.5, 
                          0.5*(c(param$gamma1[o],inits[[ch]]$gamma[[k]][[o]][1:(ncat[o]-3)]) + inits[[ch]]$gamma[[k]][[o]]), 
                          inits[[ch]]$gamma[[k]][[o]][ncat[o]-2] + 0.5)
              }else{
                Vals <- c(param$gamma1[o] - 0.5, 
                          0.5*(param$gamma1[o] + inits[[ch]]$gamma[[k]][[o]]), 
                          inits[[ch]]$gamma[[k]][[o]][ncat[o]-2] + 0.5)
              }
              inits[[ch]]$latent[[o]][selection[[k]]] <- Vals[Ordasnum[[o]][selection[[k]]]]
            }
          }else{
            probs <- cumsum(table(Y[,o]))/N
            tildemu <- param$gamma1[o]-qnorm(probs[1])
            inits[[ch]]$gamma[[o]] <- tildemu + qnorm(probs[2:(ncat[o]-1)])
            if(ncat[o]>3){
              Vals <- c(param$gamma1[o] - 0.5, 
                        0.5*(c(param$gamma1[o],inits[[ch]]$gamma[[o]][1:(ncat[o]-3)]) + inits[[ch]]$gamma[[o]]), 
                        inits[[ch]]$gamma[[o]][ncat[o]-2] + 0.5)
            }else{
              Vals <- c(param$gamma1[o] - 0.5, 
                        0.5*(param$gamma1[o] + inits[[ch]]$gamma[[o]]), 
                        inits[[ch]]$gamma[[o]][ncat[o]-2] + 0.5)
            }
            inits[[ch]]$latent[[o]] <- Vals[Ordasnum[[o]]]
          }
        } # end of else of if InitType == 3
      } # end of for o in Ords
      for(b in Bins){
        inits[[ch]]$latent[[b]] <- Y[,b]+param$gammaBin-0.5
      }
      # fixed betas 
      inits[[ch]]$beta <- list()
      if(spec["beta"]){for(k in 1:K){inits[[ch]]$beta[[k]] <- list()}}
      for(y in c(Nums,Ords,Bins)){
        Xvars <- Formula[[y]]$fixed
        if(spec["beta"]){
          for(k in 1:K){
            inits[[ch]]$beta[[k]][[y]] <- numeric(nfix[y])
            if(InitType[ch] == 1){
              if(is.element(y,Nums)){pomY <- Y[selection[[k]],y]
              }else{pomY <- inits[[ch]]$latent[[y]][selection[[k]]]}
              boldX <- as.matrix(X[selection[[k]],unique(c("intercept",Xvars))])
              coefs <- as.numeric(solve(t(boldX) %*% boldX, t(boldX) %*% pomY))
              if(is.element("intercept", Xvars)){
                inits[[ch]]$beta[[k]][[y]][Xvars == "intercept"] <- coefs[1]
                inits[[ch]]$beta[[k]][[y]][Xvars != "intercept"] <- coefs[-1]
              }else{
                inits[[ch]]$beta[[k]][[y]] <- coefs[-1]
              }
            }else{ # it is zero or randomly distributed around zero
              if(InitType[ch] == 3){inits[[ch]]$beta[[k]][[y]] <- rnorm(nfix[y], 0, param$sd_beta)}
            }
            names(inits[[ch]]$beta[[k]][[y]]) <- Xvars
          } # end of for(k in 1:K)
        }else{
          inits[[ch]]$beta[[y]] <- numeric(nfix[y])
          if(InitType[ch] == 1){
            if(is.element(y,Nums)){pomY <- Y[,y]
            }else{pomY <- inits[[ch]]$latent[[y]]}
            boldX <- as.matrix(X[,unique(c("intercept",Xvars))])
            coefs <- as.numeric(solve(t(boldX) %*% boldX, t(boldX) %*% pomY))
            if(is.element("intercept", Xvars)){
              inits[[ch]]$beta[[y]][Xvars == "intercept"] <- coefs[1]
              inits[[ch]]$beta[[y]][Xvars != "intercept"] <- coefs[-1]
            }else{
              inits[[ch]]$beta[[y]] <- coefs[-1]
            }
          }else{ # it is zero or randomly distributed around zero
            if(InitType[ch] == 3){inits[[ch]]$beta[[y]] <- rnorm(nfix[y], 0, param$sd_beta)}
          } # end of else InitType == 1
          names(inits[[ch]]$beta[[y]]) <- Xvars
        } # end of else spec["beta"]
      } # end of for y in c(Nums,Ords,Bins)
      
      # random effects
      inits[[ch]]$b <- data.frame(pom=1)
      for(y in c(Nums,Ords,Bins)){
        pom <- as.data.frame(matrix(0,ncol=nran[y],nrow=n))
        if(InitType[ch] == 1){
          if(is.element(y,Nums)){pomY <- Y[,y]
          }else{pomY <- inits[[ch]]$latent[[y]]}
          if(spec["beta"]){
            for(k in 1:K){
              pomY[selection[[k]]] <- pomY[selection[[k]]] - as.numeric(as.matrix(X[selection[[k]],Formula[[y]]$fixed]) %*% matrix(inits[[ch]]$beta[[k]][[y]],ncol=1))
            }
          }else{pomY <- pomY - as.numeric(as.matrix(X[,Formula[[y]]$fixed]) %*% matrix(inits[[ch]]$beta[[y]],ncol=1))
          }
          for(j in 1:n){
            if(nsubj[as.character(UniqSubj[j])]>nran[y]){
              coefs <- as.numeric(solve(ZtZ[[j]][[y]], tZ[[j]][[y]] %*% pomY[Y[,Id]==UniqSubj[j]]))
              pom[j,] <- coefs
            } # else it remains to be zero
          }
        }else{ 
          if(InitType[ch] == 3){
            pom <- as.data.frame(matrix(rnorm(nran[y]*n,0,param$sd_b),ncol=nran[y],nrow=n))
          }
        }# end of else of InitType[ch] == 1
        inits[[ch]]$b <- cbind(inits[[ch]]$b, pom)
      }
      inits[[ch]]$b <- as.matrix(inits[[ch]]$b[,-1])
      colnames(inits[[ch]]$b) <- rannames
      
      # tau - only for Nums
      if(spec["tau"]){
        inits[[ch]]$tau <- list()
        for(k in 1:K){
          inits[[ch]]$tau[[k]] <- numeric(length(Nums))
          names(inits[[ch]]$tau[[k]]) <- Nums
          for(y in Nums){ # just Nums, for Ords it is 1 by default
            #if(spec["beta"]){vectorbetas <- matrix(inits[[ch]]$beta[[k]][[y]],ncol=1)
            #           }else{vectorbetas <- matrix(inits[[ch]]$beta[[y]],ncol=1)}
            # if(spec["tau"]) --> then also spec["beta"]
            pomY <- Y[selection[[k]],y] - as.numeric(as.matrix(X[selection[[k]],Formula[[y]]$fixed]) %*% matrix(inits[[ch]]$beta[[k]][[y]],ncol=1)) # vectorbetas
            for(x in Formula[[y]]$random){
              pomY <- pomY - X[selection[[k]],x]*inits[[ch]]$b[NumUniqSubj[is.element(Y[,Id], Nk[[k]])],paste(y,x,sep="_")]
            }
            inits[[ch]]$tau[[k]][y] <- sum(selection[[k]])/sum(pomY^2)
          }
        } # end of for(k in 1:K)
      }else{
        inits[[ch]]$tau <- numeric(length(Nums))
        names(inits[[ch]]$tau) <- Nums
        for(y in Nums){ # just Nums, for Ords it is 1 by default
          pomY <- numeric(N)
          if(spec["beta"]){
            for(k in 1:K){
              pomY[selection[[k]]] <- Y[selection[[k]],y] - as.numeric(as.matrix(X[selection[[k]],Formula[[y]]$fixed]) %*% matrix(inits[[ch]]$beta[[k]][[y]],ncol=1))
            }
          }else{
            pomY <- Y[,y] - as.matrix(X[,Formula[[y]]$fixed]) %*% matrix(inits[[ch]]$beta[[y]],ncol=1)
          }
          for(x in Formula[[y]]$random){
            pomY <- pomY - X[,x]*inits[[ch]]$b[NumUniqSubj[as.character(Y[,Id])],paste(y,x,sep="_")]
          }
          inits[[ch]]$tau[y] <- N/sum(pomY^2)
        }
      } # end of else of if spec["tau"]
      # mu
      if(spec["mu"]){
        inits[[ch]]$mu <- list()
        for(k in 1:K){
          inits[[ch]]$mu[[k]] <- apply(inits[[ch]]$b[NumUniqSubj[as.character(Nk[[k]])],], 2, mean)
          names(inits[[ch]]$mu[[k]]) <- rannames
        }
      }else{
        inits[[ch]]$mu <- apply(inits[[ch]]$b, 2, mean)
        names(inits[[ch]]$mu) <- rannames
      } # end of else if(spec["mu"])
      # rantau  
      inits[[ch]]$rantau <- list()
      if(spec["rantau"]){
        for(k in 1:K){
          inits[[ch]]$rantau[[k]] <- list()
          for(y in c(Nums, Ords, Bins)){
            inits[[ch]]$rantau[[k]][[y]] <- rep(param$aran / param$bran, nran[y])
            names(inits[[ch]]$rantau[[k]][[y]]) <- lrannames[[y]]
          }
        }
      }else{
        for(y in c(Nums, Ords, Bins)){
          #inits[[ch]]$rantau[[y]] <- 1/Vars[lrannames[[y]]]
          # "estimate" rantau based on its prior specification/that can be totally wrong
          inits[[ch]]$rantau[[y]] <- rep(param$aran / param$bran, nran[y])
          names(inits[[ch]]$rantau[[y]]) <- lrannames[[y]]
          #??? how to estimate prior variance of mu out of data?
        }
      } # end of else if(spec["rantau"])
      
      # InvQ and InvSigma - estimate it by covariance matrix of estimated random effects
      if(spec["InvSigma"]){
        inits[[ch]]$InvSigma <- list()
        for(k in 1:K){
          if(InitType[ch] == 2){
            inits[[ch]]$InvSigma[[k]] <- diag(totnran) / rgamma(totnran, 1/param$sd_b^2)
          }else{
            inits[[ch]]$InvSigma[[k]] <- chol2inv(chol(cov(inits[[ch]]$b[NumUniqSubj[as.character(Nk[[k]])],])))
          }
        } # end of for(k in 1:K)
      }else{
        if(InitType[ch] == 2){
          inits[[ch]]$InvSigma <- diag(totnran) / rgamma(totnran, 1/param$sd_b^2)
        }else{
          inits[[ch]]$InvSigma <- chol2inv(chol(cov(inits[[ch]]$b)))
        }
      } # end of if(spec["InvSigma"])
      
      # E InvSigma = E W(Q,nu_0) = nu_0 * Q
      # InvQ approx Inv(InvSigma/nu_0)
      if(spec["InvQ"]){
        inits[[ch]]$InvQ <- list()
        for(k in 1:K){
          if (InitType[ch]==2){
            # Automatically also spec["InvSigma"] = T
            inits[[ch]]$InvQ[[k]] <- chol2inv(chol(inits[[ch]]$InvSigma[[k]]/param$nu_0))  
          }else{
            inits[[ch]]$InvQ[[k]] <- cov(inits[[ch]]$b)*param$nu_0
          } # end of else if(InitType[ch]==2)
        } # end of for(k in 1:K)
        
      }else{
        if (InitType[ch]==2){
          if(spec["InvSigma"]){
            meanInvSigma <- matrix(0, ncol = totnran, nrow = totnran)
            for(k in 1:K){meanInvSigma <- meanInvSigma + inits[[ch]]$InvSigma[[k]]/K}
            inits[[ch]]$InvQ <- chol2inv(chol(meanInvSigma/param$nu_0))
          }else{
            inits[[ch]]$InvQ <- chol2inv(chol(inits[[ch]]$InvSigma/param$nu_0))
          } # end of if(spec["InvSigma"])
        }else{
          inits[[ch]]$InvQ <- cov(inits[[ch]]$b)*param$nu_0
        } # end of else if(InitType[ch]==2)
      } # end of else if(spec["InvQ"])
    } # end of ch in 1:Nchains
  }else{ # end of if (missing(inits))
    # Make sure that Nchains corresponds to the given number of inits
  }
  
  if(howsave == "matrix"){
    mcmc$all <- data.frame()
  }
  
  cat("Inits are ready.\n")
  
  ### Preparations for transfer from C to List output
  # creating settings matrix for all (calculable) parameters
  params <- c("w", "U", "pUik", "pUik_nonb", "latent", "tau",  
              "mu", "InvSigma", "detInvSigma", "detSigma", "Sigma", 
              "sdSigma", "corSigma", "InvQ", "detInvQ", "detQ", "Q", "b")
  otherparam <- c("gamma", "min", "max", "rantau", "beta")
  
  settings <- data.frame(save = sapply(c(params,otherparam), function(p){c(na.omit(c(whatsave[p], calc[p])))}),
                         isspec = sapply(c(params,otherparam), function(p){ifelse(is.na(spec[p]), F, spec[p])}),
                         K = K,
                         isy = F,
                         ynums = 0,
                         yords = 0,
                         ybins = 0,
                         BM = rep(B+M, length(c(params,otherparam))),
                         d1 = 0,
                         d2 = 0, 
                         BYROW = T,
                         sym = 0,
                         diag = 0,
                         diagval = 1,
                         D = 0)
  
  rownames(settings) <- c(params,otherparam)
  settings
  
  # Individual changes
  settings["w", c("d1", "D")] = c(K, 1)
  settings["U", c("d1", "D")] = c(n, 1)
  settings["pUik", c("d1", "d2", "BYROW", "D")] = c(n, K, F, 2)
  settings["pUik_nonb", c("d1", "d2", "BYROW", "D")] = c(n, K, F, 2)
  
  # otherparams (gamma, min, max) have to be done separately
  # due to d1 varying on y
  settings["gamma", c("isy", "yords", "D")] = c(T,1,1)
  settings["min", c("isy", "yords", "D")] = c(T,1,1)
  settings["max", c("isy", "yords", "D")] = c(T,1,1)
  
  settings["latent", c("isy", "yords", "ybins", "d1", "D")] = c(T,1,1,N,1)
  settings["tau", c("isy", "ynums")] = c(T, 1)
  
  # rantau and beta have d1 varying on y
  settings["rantau", c("isy", "ynums", "yords", "ybins", "D")] = c(T,1,1,1,1)
  settings["beta", c("isy", "ynums", "yords", "ybins", "D")] = c(T,1,1,1,1)
  
  settings["mu", c("d1", "D")] = c(totnran, 1)
  
  settings["InvSigma", c("d1","d2", "sym", "diag", "D")] = c(totnran,totnran,1,1,2)
  settings["detInvSigma", c("isspec")] = c(spec["InvSigma"])
  settings["detSigma", c("isspec")] = c(spec["InvSigma"])
  settings["Sigma", c("isspec","d1","d2", "sym", "diag", "D")] = c(spec["InvSigma"],totnran,totnran,1,1,2)
  settings["sdSigma", c("isspec", "d1", "D")] = c(spec["InvSigma"], totnran, 1)
  settings["corSigma", c("isspec","d1","d2", "sym", "diag", "diagval", "D")] = c(spec["InvSigma"],totnran,totnran,1,0,1,2)
  
  settings["InvQ", c("d1","d2", "sym", "diag", "D")] = c(totnran,totnran,1,1,2)
  settings["detInvQ", c("isspec")] = c(spec["InvQ"])
  settings["detQ", c("isspec")] = c(spec["InvQ"])
  settings["Q", c("isspec","d1","d2", "sym", "diag", "D")] = c(spec["InvQ"],totnran,totnran,1,1,2)
  
  settings["b", c("d1", "d2", "BYROW", "D")] = c(n, totnran, F, 2)
  yspecd1 <- list()
  yspecd1$gamma <- yspecd1$max <- yspecd1$min <- ncat-2
  yspecd1$rantau <- nran
  yspecd1$beta <- nfix
    
  ### Preparation of parameters for Gibbs sampler in C
  ## Data
  # First column is going to be the id variable (0-th column)
  cId <- NumUniqSubj[as.character(Y[,Id])] - 1
  # -1 is there for C which works better with number beggining with 0
  # other columns  (beggining with 1st column)
  cY <- cX <- numeric()
  for(y in c(Nums, Ords, Bins)){cY <- c(cY, as.numeric(as.character(Y[,y])))}
  for(x in setdiff(colnames(X), Id)){cX <- c(cX, X[,x])} # no need for Id
  # Formula
  cFormulaF <- cFormulaR <- numeric()
  for(y in c(Nums, Ords, Bins)){
    cFormulaF <- c(cFormulaF, which(is.element(setdiff(colnames(X), Id), Formula[[y]]$fixed)))
    cFormulaR <- c(cFormulaR, which(is.element(setdiff(colnames(X), Id), Formula[[y]]$random)))
  }
  cFormulaF <- cFormulaF - 1
  cFormulaR <- cFormulaR - 1
  cnfix <- as.numeric(nfix) # number of FIXED  regressors with variables y
  cnran <- as.numeric(nran) # number of RANDOM regressors with variables y
  # dims - in the following order:
  # w, pUik, U can be determined out of K and n
  cdims <- c(K, # w
             n, # U
             n*K, # pUik
             sum(ncat-2), # gamma
             sum(ncat-2), # min
             sum(ncat-2), # max
             (nY["Ords"]+nY["Bins"])*N, # latent
             totnfix, # beta
             nY["Nums"], # tau
             totnran, # mu
             totnran, # rantau
             totnran*(totnran+1)/2, # InvSigma
             totnran*(totnran+1)/2, # InvQ
             n*totnran, # b
             # calculable parameters
             n*K, # pUik_nonb
             totnran*(totnran+1)/2, # Sigma
             totnran*(totnran+1)/2, # Q
             1,1,1,1, # determinants of Sigma, InvSigma, Q, InvQ
             totnran, # sdSigma
             totnran*(totnran-1)/2 # corSigma
  )
  cdimswithK <- c(K, # w
                  n, # U
                  n*K, # pUik
                  sum(ncat-2)*ifelse(spec["gamma"], K, 1), # gamma
                  sum(ncat-2)*ifelse(spec["gamma"], K, 1), # min
                  sum(ncat-2)*ifelse(spec["gamma"], K, 1), # max
                  (nY["Ords"]+nY["Bins"])*N, # latent
                  totnfix*ifelse(spec["beta"], K, 1), # beta
                  nY["Nums"]*ifelse(spec["tau"], K, 1), # tau
                  totnran*ifelse(spec["mu"], K, 1), # mu
                  totnran*ifelse(spec["rantau"], K, 1), # rantau
                  totnran*(totnran+1)/2*ifelse(spec["InvSigma"], K, 1), # InvSigma
                  totnran*(totnran+1)/2*ifelse(spec["InvQ"], K, 1), # InvQ
                  n*totnran, # b
                  # calculable parameters
                  n*K, # pUik_nonb
                  totnran*(totnran+1)/2*ifelse(spec["InvSigma"], K, 1), # Sigma
                  totnran*(totnran+1)/2*ifelse(spec["InvQ"], K, 1), # Q
                  ifelse(spec["InvSigma"], K, 1), # detSigma
                  ifelse(spec["InvSigma"], K, 1), # detInvSigma
                  ifelse(spec["InvQ"], K, 1), # detQ
                  ifelse(spec["InvQ"], K, 1), # detInvQ
                  totnran*ifelse(spec["InvSigma"], K, 1), # sdSigma
                  totnran*(totnran-1)/2*ifelse(spec["InvSigma"], K, 1) # corSigma
  )
  names(cdims) <- names(cdimswithK) <- c("w","U","pUik","gamma","min","max","latent",
                                         "beta","tau","mu","rantau","InvSigma","InvQ","b",
                                         "pUik_nonb","Sigma","Q","detSigma","detInvSigma",
                                         "detQ","detInvQ","sdSigma","corSigma")
  for(p in names(cdims)){
    settings[p, "dims"] = cdims[p]
    settings[p, "dimswithK"] = cdimswithK[p]
  }
  # those dimensions are calculable in C function as well, when
  # K, n, ncat, nY, spec, nran, nfix     are available
  
  # cparam
  cparam <- param
  cparam$sd_beta <- cparam$sd_b <- cparam$V <- NULL
  cparam$fixD <- unlist(cparam$fixD)
  cparam$fixmu <- unlist(cparam$fixmu)
  cparam$InvV <- param$InvV[upper.tri(param$InvV, diag = TRUE)]
  # order is nu_0, nu_1, delta, a, b, aran, bran, (all [1])
  # ranmu0 [totnran], 
  # gamma1 [nY[Ord]], 
  # gammaBin [1], 
  # fixD [totnfix],
  # fixmu [totnran],
  # InvV [totnran*(totnran+1)/2]
  cparam <- unlist(cparam) # will be delivered as vector
  
  # arrays/fields where those variables will (might) be stored
    if(whatsave[1]){cw          <- double((B+M)*cdimswithK[1])}else{cw <- as.double(0)}
    if(whatsave[2]){cU          <- integer((B+M)*cdimswithK[2])}else{cU <- as.integer(0)}
    if(whatsave[3]){cpUik       <- double((B+M)*cdimswithK[3])}else{cpUik <- as.double(0)}
    if(whatsave[4]){cgamma      <- double((B+M)*cdimswithK[4])}else{cgamma <- as.double(0)}
    if(whatsave[5]){cmin        <- double((B+M)*cdimswithK[5])}else{cmin <- as.double(0)}
    if(whatsave[6]){cmax        <- double((B+M)*cdimswithK[6])}else{cmax <- as.double(0)}
    if(whatsave[7]){clatent     <- double((B+M)*cdimswithK[7])}else{clatent <- as.double(0)}
    if(whatsave[8]){cbeta       <- double((B+M)*cdimswithK[8])}else{cbeta <- as.double(0)}
    if(whatsave[9]){ctau        <- double((B+M)*cdimswithK[9])}else{ctau <- as.double(0)}
    if(whatsave[10]){cmu        <- double((B+M)*cdimswithK[10])}else{cmu <- as.double(0)}
    if(whatsave[11]){crantau    <- double((B+M)*cdimswithK[11])}else{crantau <- as.double(0)}
    if(whatsave[12]){cInvSigma  <- double((B+M)*cdimswithK[12])}else{cInvSigma <- as.double(0)}
    if(whatsave[13]){cInvQ      <- double((B+M)*cdimswithK[13])}else{cInvQ <- as.double(0)}
    if(whatsave[14]){cb         <- double((B+M)*cdimswithK[14])}else{cb <- as.double(0)}
    if(calc[1]){     cpUik_nonb <- double((B+M)*cdimswithK[15])}else{cpUik_nonb <- as.double(0)}
    if(calc[2]){     cQ         <- double((B+M)*cdimswithK[17])}else{cQ <- as.double(0)}
    if(calc[3]){     cdetQ      <- double((B+M)*cdimswithK[20])}else{cdetQ <- as.double(0)}
    if(calc[4]){     cdetInvQ   <- double((B+M)*cdimswithK[21])}else{cdetInvQ <- as.double(0)}
    if(calc[5]){     cSigma     <- double((B+M)*cdimswithK[16])}else{cSigma <- as.double(0)}
    if(calc[6]){     cdetSigma  <- double((B+M)*cdimswithK[18])}else{cdetSigma <- as.double(0)}
    if(calc[7]){  cdetInvSigma  <- double((B+M)*cdimswithK[19])}else{cdetInvSigma <- as.double(0)}
    if(calc[8]){     csdSigma   <- double((B+M)*cdimswithK[22])}else{csdSigma<- as.double(0)}
    if(calc[9]){     ccorSigma  <- double((B+M)*cdimswithK[23])}else{ccorSigma<- as.double(0)}
    
  # dimension of all ZtZ
    cdimZtZ <- sum(cnran * (cnran+1)/2) 
    
  # list of last generated states
  last <- list()
    
  ### Now prepare inits for each chain
  # and generate chains
  for(ch in 1:Nchains){
    initsch <- inits[[ch]]
    # cinits
    cinits <- list()
    cinits$w <- initsch$w
    # cinits$U <- initsch$U # those are integers will be passed separately
    cinits$pUik <- c(initsch$pUik)
    cinits$gamma <- unlist(initsch$gamma)
    cinits$latent <- unlist(initsch$latent)
    cinits$beta <- unlist(initsch$beta)
    cinits$tau <- unlist(initsch$tau)
    cinits$mu <- unlist(initsch$mu)
    cinits$rantau <- unlist(initsch$rantau)
    if(spec["InvSigma"]){
      InvSigma2 <- list()
      for(k in 1:K){
        InvSigma2[[k]] <- initsch$InvSigma[[k]][upper.tri(initsch$InvSigma[[k]], diag = T)]
      }
      cinits$InvSigma <- unlist(InvSigma2)
    }else{
      cinits$InvSigma <- initsch$InvSigma[upper.tri(initsch$InvSigma, diag = T)]
    }
    if(spec["InvQ"]){
      InvQ2 <- list()
      for(k in 1:K){
        InvQ2[[k]] <- initsch$InvQ[[k]][upper.tri(initsch$InvQ[[k]], diag = T)]
      }
      cinits$InvQ <- unlist(InvQ2)
    }else{
      cinits$InvQ <- initsch$InvQ[upper.tri(initsch$InvQ, diag = T)]
    }
    cinits$b <- c(initsch$b)
    
    cinits <- unlist(cinits) # will be delivered as a vector
    cinits_U <- initsch$U # expects values 0, ..., K-1
    
    #dyn.load(paste(gsub("EU-SILC", "Classification", WD), "class_num_ord_bin.dll", sep = "/"))
    #dyn.unload(paste(gsub("EU-SILC", "Classification", WD), "class_num_ord_bin.dll", sep = "/"))
    
    
    #system.time(
      cmcmc <-
        .C("gibbs_class_num_ord_bin",
           Id        = as.integer(cId),
           Y         = as.double(cY),
           X         = as.double(cX),
           spec      = as.integer(spec),
           calc      = as.integer(calc),
           whatsave  = as.integer(whatsave),
           param     = as.double(cparam), # passed as vector of values 
           vecinits  = as.double(cinits), # passed as vector of double values (except U)
           inits_U   = as.integer(cinits_U), # passed as vector of integers
           veclast   = double(length(cinits)), # passed as vector of double values (except U)
           last_U    = integer(length(cinits_U)), # passed as vector of integers
           # parameters describing dimensions
           chain     = as.integer(ch), # number of the chain
           K         = as.integer(K), # number of classes
           BM        = as.integer(B + M), # total number of generated states
           FormulaF  = as.integer(cFormulaF), # numbers of columns of X that should be used for FIXED  effects of modelled responses
           FormulaR  = as.integer(cFormulaR), # numbers of columns of X that should be used for RANDOM effects of modelled responses
           N         = as.integer(N), # total number of observations
           dims      = as.integer(cdims), # the length of subarray that corresponds to one state (disected by various parameters)
           dimswithK = as.integer(cdimswithK), # the length of subarray that corresponds to one state 
           # (disected by various parameters, also multiplication by K incorporated when such parameters is class-specific)
           ncat      = as.integer(ncat), # the counts of categories of ordinal variables
           nY        = as.integer(nY), # 3 numbers: counts of Nums, Ords and Bins variables
           n         = as.integer(n), # total number of subjects (different ids in the dataset)
           nfix      = as.integer(cnfix), # 
           nran      = as.integer(cnran),
           dimZtZ    = as.integer(cdimZtZ),
           #ZtZ       = double(n*cdimZtZ),
           #predictor = double(N*sum(nY)),
           # the function should count totnran, totnfix and cumsum versions of nfix and nran
           # arrays to store generated states
           w           = cw,
           U           = cU,
           pUik        = cpUik,
           latent      = clatent,
           gamma       = cgamma,
           min         = cmin,
           max         = cmax,
           beta        = cbeta,
           tau         = ctau,
           mu          = cmu,
           rantau      = crantau,
           InvSigma    = cInvSigma,
           InvQ        = cInvQ,
           b           = cb,
           pUik_nonb   = cpUik_nonb,
           Sigma       = cSigma,
           Q           = cQ,
           detSigma    = cdetSigma,
           detInvSigma = cdetInvSigma,
           detQ        = cdetQ,
           detInvQ     = cdetInvQ,
           sdSigma     = csdSigma,
           corSigma    = ccorSigma
        )
    #)
    cat(paste0("Sampling of chain ", ch, " is completed.\n"))
      
    ### reconstruction of last state from cmcmc$last_U and cmcmc$veclast
    last[[ch]] <- list()
    # order of construction of inits: w, U, pUik, gamma, latent, beta, 
                                      # b, tau, mu, rantau, InvSigma, InvQ
    # order of parameters in veclast: w, pUik, gamma, latent, beta, tau, mu, 
                                      # rantau, InvSigma, InvQ, b
    if(howsave != "cmcmc"){
      last[[ch]]$w <- cmcmc$veclast[1:(lastcopied <- cdims[1])]
      last[[ch]]$U <- cmcmc$last_U
      last[[ch]]$pUik <- matrix(cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdimswithK[3])],
                                n, K, byrow = F)
      # gamma
      if(spec["gamma"]){
        last[[ch]]$gamma <- list()
        for(k in 1:K){
          last[[ch]]$gamma[[k]] <- list()
          for(o in Ords){
            last[[ch]]$gamma[[k]][[o]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + ncat[o]-2)]
          }
        }
      }else{
        last[[ch]]$gamma <- list()
        for(o in Ords){
          last[[ch]]$gamma[[o]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + ncat[o]-2)]
        }
      }
      # latent
      last[[ch]]$latent <- list()
      for(y in c(Ords, Bins)){
        last[[ch]]$latent[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + N)]
      }
      # beta
      if(spec["beta"]){
        last[[ch]]$beta <- list()
        for(k in 1:K){
          last[[ch]]$beta[[k]] <- list()
          for(y in c(Nums, Ords, Bins)){
            last[[ch]]$beta[[k]][[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nfix[y])]
          }
        }
      }else{
        last[[ch]]$beta <- list()
        for(y in c(Nums, Ords, Bins)){
          last[[ch]]$beta[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nfix[y])]
        }
      }
      # tau
      if(spec["tau"]){
        last[[ch]]$tau <- list()
        for(k in 1:K){
          last[[ch]]$tau[[k]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[9])]
        }
      }else{
        last[[ch]]$tau <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[9])]
      }
      # mu
      if(spec["mu"]){
        last[[ch]]$mu <- list()
        for(k in 1:K){
          last[[ch]]$mu[[k]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[10])]
        }
      }else{
        last[[ch]]$mu <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[10])]
      }
      # rantau
      if(spec["rantau"]){
        last[[ch]]$rantau <- list()
        for(k in 1:K){
          last[[ch]]$rantau[[k]] <- list()
          for(y in c(Nums, Ords, Bins)){
            last[[ch]]$rantau[[k]][[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nran[y])]
          }
        }
      }else{
        last[[ch]]$rantau <- list()
        for(y in c(Nums, Ords, Bins)){
          last[[ch]]$rantau[[y]] <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + nran[y])]
        }
      }
      # InvSigma
      if(spec["InvSigma"]){
        last[[ch]]$InvSigma <- list()
        for(k in 1:K){
          pommatrix <- matrix(0, totnran, totnran)
          pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[12])]
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          pommatrix <- t(pommatrix)
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          last[[ch]]$InvSigma[[k]] <- pommatrix
        }
      }else{
        pommatrix <- matrix(0, totnran, totnran)
        pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[12])]
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        pommatrix <- t(pommatrix)
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        last[[ch]]$InvSigma <- pommatrix
      }
      # InvQ
      if(spec["InvQ"]){
        last[[ch]]$InvQ <- list()
        for(k in 1:K){
          pommatrix <- matrix(0, totnran, totnran)
          pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[13])]
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          pommatrix <- t(pommatrix)
          pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
          last[[ch]]$InvQ[[k]] <- pommatrix
        }
      }else{
        pommatrix <- matrix(0, totnran, totnran)
        pomvec <- cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdims[13])]
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        pommatrix <- t(pommatrix)
        pommatrix[upper.tri(pommatrix, diag = T)] <- pomvec
        last[[ch]]$InvQ <- pommatrix
      }
      # b
      last[[ch]]$b <- matrix(cmcmc$veclast[(lastcopied + 1):(lastcopied <- lastcopied + cdimswithK[14])],
                             n, totnran, byrow = F)
    }
    cat(paste0("Saving the last state of chain ", ch, " is completed.\n"))
    
    ### Saving results as structured list (converted from cmcmc)
    # uses previously created matrix settings
    if(howsave == "cmcmc"){
      mcmc[[ch]] <- cmcmc
    }
    
    if(howsave == "list"){
      # results will be returned in structured list - the same way as original R function
      chain <- list()
      for(p in c(params,otherparam)){
        if(settings[p, "save"]){
          if(is.element(p, otherparam)){
            chain[[p]] <- FromCtoList_settings(values = cmcmc[[p]], 
                                               p = p, 
                                               settings = settings,
                                               yspecd1 = yspecd1[[p]])
          }else{
            chain[[p]] <- FromCtoList_settings(values = cmcmc[[p]], 
                                               p = p, 
                                               settings = settings)
          }
        }
      }
      
      mcmc[[ch]] <- chain
      mcmc[[ch]]$inits <- inits[[ch]]
      mcmc[[ch]]$last <- last[[ch]]
      if(!initsgiven){
        mcmc[[ch]]$InitType <- InitType[ch]
      }else{
        mcmc[[ch]]$InitType <- "given"
      }
      cat(paste0("Saving the chain ", ch, " into lists is completed.\n"))
    }
    
    if(howsave == "matrix"){
      # results will be returned in matrix (row = state number, col = variable)
      # column chain will distinguish from what chains it comes from
      
      AllData <- matrix(ch, nrow = B+M, ncol = 1)
      colnames(AllData) <- "chain"
      
      for(p in c(params,otherparam)){
        if(settings[p, "save"]){
          if(is.element(p, otherparam)){
            AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[p]],
                                                             p = p,
                                                             settings = settings,
                                                             yspecd1 = yspecd1[[p]],
                                                             Nums = Nums, Ords = Ords, Bins = Bins))
          }else{
            AllData <- cbind(AllData, FromCtoMatrix_settings(values = cmcmc[[p]],
                                                             p = p,
                                                             settings = settings,
                                                             Nums = Nums, Ords = Ords, Bins = Bins))
          }
        }
      }
      
      mcmc$all <- rbind(mcmc$all, AllData)
      mcmc$inits <- inits
      if(!initsgiven){
        mcmc$InitType <- InitType
      }else{
        mcmc$InitType <- "given"
      }
      cat(paste0("Saving the chain ", ch, " into one big matrix is completed.\n"))
    }
    
  } # end of chain
  cat("All chains have been already sampled.\n")
  
  mcmc$M <- M
  mcmc$B <- B
  mcmc$Nchains <- Nchains
  mcmc$BM <- B + M
  mcmc$Nums <- Nums
  mcmc$Ords <- Ords
  mcmc$Bins <- Bins
  mcmc$Formula <- Formula
  mcmc$K <- K
  mcmc$spec <- spec
  mcmc$calc <- calc
  mcmc$whatsave <- whatsave
  mcmc$ncat <- ncat
  mcmc$howsave <- howsave
  mcmc$last <- last
  mcmc$param <- param
  mcmc$nfix <- nfix
  mcmc$nran <- nran
  mcmc$totnfix <- totnfix
  mcmc$totnran <- totnran
  mcmc$lfixnames <- lfixnames
  mcmc$lrannames <- lrannames
  mcmc$fixnames <- fixnames
  mcmc$rannames <- rannames
  mcmc$nY <- nY
  mcmc$settings <- settings
  mcmc$yspecd1 <- yspecd1
  return(mcmc)
} # end of function