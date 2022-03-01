### All plotting functions
### - descriptive for generated dataset
### - traceplot, density, ACF, ECDF
### - classes comparison
library("colorspace")

PlotPanelDataIntoOne <- function(data, yvar, K, addjitter = F,
                                 subjects = 1:n){
  col <- rainbow_hcl(K, c = 80, l =50)
  collight <- rainbow_hcl(K, c = 80, l =70)
  
  if(addjitter){
    data$resp <- as.numeric(as.character(data[[yvar]])) + runif(dim(data)[1], -0.3, 0.3)
  }else{
    data$resp <- data[[yvar]]
  }
  plot(x = c(0,1), y = c(0,1), type = "n", 
       xlim = c(0, 1),
       ylim = c(min(data$resp), max(data$resp)),
       xlab = "Time",
       ylab = yvar)
  for (i in subjects){
    datai <- subset(data, subject == i)
    lines(datai$time, datai$resp, col = collight[datai$U[1]])
  }
  for (k in 1:K){
    lmmodel <- lm(resp ~ time, data[data$U==k,])
    abline(lmmodel, col = col[k], lwd = 2)
  }
  legend("topleft", legend = 1:K, title = "Class", lty = 1, lwd = 2,
         col = col, pt.bg = collight, bty = "n")
}

PlotPanelDataClassified <- function(data, yvar, K, addjitter = F,
                                    subjects = 1:n, class){
  col <- c("grey40", rainbow_hcl(K, c = 80, l =50))
  collight <- c("grey70", rainbow_hcl(K, c = 80, l =70))
  
  if(addjitter){
    data$resp <- as.numeric(as.character(data[[yvar]])) + runif(dim(data)[1], -0.3, 0.3)
  }else{
    data$resp <- data[[yvar]]
  }
  plot(x = c(0,1), y = c(0,1), type = "n", 
       xlim = c(0, 1),
       ylim = c(min(data$resp), max(data$resp)+ifelse(addjitter,1,3)),
       xlab = "Time",
       ylab = yvar)
  for (i in subjects){
    datai <- subset(data, subject == i)
    lines(datai$time, datai$resp, 
          col = ifelse(class[i] == datai$U[1], collight[class[i]+1], col[class[i]+1]), 
          lty = ifelse(class[i] == datai$U[1], 1, 2),
          lwd = ifelse(class[i] == datai$U[1], 1, 2))
  }
  
  legend("topleft", legend = c("None", 1:K), title = "Class", 
         lty = c(2, rep(1,K)), lwd = 2,
         col = c(col[1], collight[1:K + 1]),  bty = "n")
  legend("topright", 
         legend = c(paste0("Matches the true class (", 
                           format(100*sum(class == data$U[seq(1,dim(data)[1],by = dim(datai)[1])])/n, 
                                  digits = 2, nsmall = 2), " %)"), 
                    paste0("Un/misclassified (", 
                           format(100*sum(class != data$U[seq(1,dim(data)[1],by = dim(datai)[1])])/n, 
                                  digits = 2, nsmall = 2), " %)")), 
         title = "Correctness", lty = c(1,2), lwd = c(1,2),
         col = c(collight[1], col[1]),  bty = "n")
}


### Computing posterior estimates
library("coda")
# Creating list of mcmc tables 
# each row represents one iteration of Gibbs
# each column represents one of the variables
# selected - which variables are about to be used 
#          - uses general names -> includes the whole block of variables
CreateAllData <- function(mcmc, selected, B, M, thin, ChainsToUse){
  if (missing(B)){B <- mcmc$B}
  if (missing(M)){M <- mcmc$BM - B}
  if (missing(thin)){thin <- 1}
  if (missing(ChainsToUse)){ChainsToUse <- 1:mcmc$Nchains}
  AllData <- list()
  for(ch in ChainsToUse){
    datach <- data.frame(pom=1)
    # parameters concerning Classification
    if(is.element("U", selected)){
      addtodatach <- as.data.frame(mcmc[[ch]]$U[seq(B+1,B+M,by=thin),])
      colnames(addtodatach) <- paste0("U[",1:length(unique(Y[,Id])),"]")
      datach <- cbind(datach, addtodatach)
    }
    if(is.element("pUik", selected)){
      for(k in 1:mcmc$K){
        addtodatach <- as.data.frame(mcmc[[ch]]$pUik[seq(B+1,B+M,by=thin),,k])
        colnames(addtodatach) <- paste0("pUik_",k,"[",1:length(unique(Y[,Id])),"]")
        datach <- cbind(datach, addtodatach)
      }
    }
    if(is.element("w", selected)){
      addtodatach <- as.data.frame(mcmc[[ch]]$w[seq(B+1,B+M,by=thin),])
      colnames(addtodatach) <- paste0("w[",1:mcmc$K,"]")
      datach <- cbind(datach, addtodatach)
    }
    # parameters concerning Ords and latent parametrization
    for(ovar in c("gamma", "min", "max")){
      if(is.element(ovar, selected)){
        if(mcmc$spec["gamma"]){
          for(o in mcmc$Ords){
            for(k in 1:mcmc$K){
              addtodatach <- as.data.frame(mcmc[[ch]][[ovar]][[k]][[o]][seq(B+1,B+M,by=thin),])
              colnames(addtodatach) <- paste0(ovar,"_",o,"(",k,")[",1:(ncol(mcmc[[ch]][[ovar]][[k]][[o]])),"]")
              datach <- cbind(datach, addtodatach)
            }
          }
        }else{
          for(o in mcmc$Ords){
            addtodatach <- as.data.frame(mcmc[[ch]][[ovar]][[o]][seq(B+1,B+M,by=thin),])
            colnames(addtodatach) <- paste0(ovar,"_",o,"[",1:(ncol(mcmc[[ch]][[ovar]][[o]])),"]")
            datach <- cbind(datach, addtodatach)
          }
        }
      }
    } # end of for ovar in gamma, min, max
    if(is.element("latent", selected)){
      for(o in mcmc$Ords){
        addtodatach <- as.data.frame(mcmc[[ch]]$latent[[o]][seq(B+1,B+M,by=thin),])
        colnames(addtodatach) <- paste0(ovar,"_",o,"[",1:(ncol(mcmc[[ch]]$latent[[o]])),"]")
        datach <- cbind(datach, addtodatach)
      }
    }
    # tau
    if(is.element("tau",selected)){
      if(mcmc$spec["tau"]){
        for(k in 1:mcmc$K){
          addtodatach <- as.data.frame(mcmc[[ch]]$tau[[k]][seq(B+1,B+M,by=thin),])
          colnames(addtodatach) <- paste0("tau_",mcmc$Nums,"(",k,")")
          datach <- cbind(datach, addtodatach)
        }
      }else{
        addtodatach <- as.data.frame(mcmc[[ch]]$tau[seq(B+1,B+M,by=thin),])
        colnames(addtodatach) <- paste0("tau_",mcmc$Nums)
        datach <- cbind(datach, addtodatach)
      }
    }
    # rantau
    if(is.element("rantau",selected)){
      if(mcmc$spec["rantau"]){
        for(k in 1:mcmc$K){
          for(y in c(mcmc$Nums,mcmc$Ords,mcmc$Bins)){
            addtodatach <- as.data.frame(mcmc[[ch]]$rantau[[k]][[y]][seq(B+1,B+M,by=thin),])
            colnames(addtodatach) <- paste0("rantau_",y,"(",k,")[",mcmc$formula[[y]]$random,"]")
            datach <- cbind(datach, addtodatach)
          }
        }
      }else{
        for(y in c(mcmc$Nums,mcmc$Ords,mcmc$Bins)){
          addtodatach <- as.data.frame(mcmc[[ch]]$rantau[[y]][seq(B+1,B+M,by=thin),])
          colnames(addtodatach) <- paste0("rantau_",y,"[",mcmc$formula[[y]]$random,"]")
          datach <- cbind(datach, addtodatach)
        }
      }
    }
    # beta
    if(is.element("beta",selected)){
      if(mcmc$spec["beta"]){
        for(k in 1:mcmc$K){
          for(y in c(Nums,Ords)){
            addtodatach <- as.data.frame(mcmc[[ch]]$beta[[k]][[y]][seq(B+1,B+M,by=thin),])
            colnames(addtodatach) <- paste0("beta_",y,"(",k,")[",mcmc$formula[[y]]$fixed,"]")
            datach <- cbind(datach, addtodatach)
          }
        }
      }else{
        for(y in c(Nums,Ords)){
          addtodatach <- as.data.frame(mcmc[[ch]]$beta[[y]][seq(B+1,B+M,by=thin),])
          colnames(addtodatach) <- paste0("beta_",y,"[",mcmc$formula[[y]]$fixed,"]")
          datach <- cbind(datach, addtodatach)
        }
      }
    }
    # rannames for matrices coordinates
    rannames <- numeric()
    for(y in c(mcmc$Nums, mcmc$Ords, mcmc$Bins)){
      rannames <- c(rannames, paste(y, Formula[[y]]$random, sep = "_"))
    }
    # mu 
    if(is.element("mu",selected)){
      if(mcmc$spec["mu"]){
        for(k in 1:mcmc$K){
          addtodatach <- as.data.frame(mcmc[[ch]]$mu[[k]][seq(B+1,B+M,by=thin),])
          colnames(addtodatach) <- paste0("mu_",rannames,"(",k,")")
          datach <- cbind(datach, addtodatach)
        }
      }else{
        addtodatach <- as.data.frame(mcmc[[ch]]$mu[seq(B+1,B+M,by=thin),])
        colnames(addtodatach) <- paste0("mu_",rannames)
        datach <- cbind(datach, addtodatach)
      }
    }
    # InvQ, InvSigma
    for(matvar in c("InvQ","InvSigma")){
      if(is.element(matvar, selected)){
        if(mcmc$spec[matvar]){
          for(k in 1:mcmc$K){
            for(i in 1:dim(mcmc[[ch]][[matvar]][[k]])[2]){
              for(j in (i+1):dim(mcmc[[ch]][[matvar]][[k]])[2]){
                addtodatach <- as.data.frame(mcmc[[ch]][[matvar]][[k]][seq(B+1,B+M,by=thin),i,j])
                colnames(addtodatach) <- paste0(matvar,"(",k,")[",rannames[i],",",rannames[j],"]")
                datach <- cbind(datach, addtodatach)
              }
            }
          }
        }else{
          for(i in 1:dim(mcmc[[ch]][[matvar]])[2]){
            for(j in (i+1):dim(mcmc[[ch]][[matvar]])[2]){
              addtodatach <- as.data.frame(mcmc[[ch]][[matvar]][seq(B+1,B+M,by=thin),i,j])
              colnames(addtodatach) <- paste0(matvar,"[",rannames[i],",",rannames[j],"]")
              datach <- cbind(datach, addtodatach)
            }
          }
        }
      }
    } # end of for matvar
    # Sigma, detSigma, detInvQ, corSigma, sdSigma
    if(is.element("Sigma", selected)){
      if(mcmc$spec["InvSigma"]){
        for(k in 1:mcmc$K){
          Sigmas <- corSigmas <- array(0, dim = c(length(seq(B+1, B+M, by = thin)), 
                                                  dim(mcmc[[ch]][["InvSigma"]][[k]])[2],
                                                  dim(mcmc[[ch]][["InvSigma"]][[k]])[2]))
          sdSigmas <- matrix(0, length(seq(B+1, B+M, by = thin)), dim(mcmc[[ch]][["InvSigma"]][[k]])[2])
          detSigma <-  numeric(length(seq(B+1, B+M, by = thin)))
          for(i in seq(B+1, B+M, by = thin)){
            Sigmas[1+(i-B-1)/thin,,] <- chol2inv(chol(mcmc[[ch]]$InvSigma[[k]][i,,]))
            sdSigmas[1+(i-B-1)/thin,] <- sqrt(diag(Sigmas[1+(i-B-1)/thin,,]))
            corSigmas[1+(i-B-1)/thin,,] <- diag(1/sdSigmas[1+(i-B-1)/thin,]) %*% Sigmas[1+(i-B-1)/thin,,] %*% diag(1/sdSigmas[1+(i-B-1)/thin,])
            detSigma[1+(i-B-1)/thin] <- det(Sigmas[1+(i-B-1)/thin,,])
          }
          for(i in 1:dim(mcmc[[ch]][["InvSigma"]][[k]])[2]){
            for(j in i:dim(mcmc[[ch]][["InvSigma"]][[k]])[2]){
              addtodatach <- as.data.frame(Sigmas[,i,j])
              colnames(addtodatach) <- paste0("Sigma(",k,")[",rannames[i],",",rannames[j],"]")
              datach <- cbind(datach, addtodatach)
            }
          }
          for(i in 1:(dim(mcmc[[ch]][["InvSigma"]][[k]])[2]-1)){
            for(j in (i+1):dim(mcmc[[ch]][["InvSigma"]][[k]])[2]){
              addtodatach <- as.data.frame(corSigmas[,i,j])
              colnames(addtodatach) <- paste0("corSigma(",k,")[",rannames[i],",",rannames[j],"]")
              datach <- cbind(datach, addtodatach)
            }
          }
          for(i in 1:dim(mcmc[[ch]][["InvSigma"]][[k]])[2]){
            addtodatach <- as.data.frame(sdSigmas[,i])
            colnames(addtodatach) <- paste0("sdSigma(",k,")[",rannames[i],"]")
            datach <- cbind(datach, addtodatach)
          }
          addtodatach <- as.data.frame(detSigma)
          colnames(addtodatach) <- paste0("detSigma(",k,")")
          datach <- cbind(datach, addtodatach)
        }
      }else{
        Sigmas <- corSigmas <- array(0, dim = c(length(seq(B+1, B+M, by = thin)), 
                                                dim(mcmc[[ch]][["InvSigma"]])[2],
                                                dim(mcmc[[ch]][["InvSigma"]])[2]))
        sdSigmas <- matrix(0, length(seq(B+1, B+M, by = thin)), dim(mcmc[[ch]][["InvSigma"]])[2])
        detSigma <-  numeric(length(seq(B+1, B+M, by = thin)))
        for(i in seq(B+1, B+M, by = thin)){
          Sigmas[1+(i-B-1)/thin,,] <- chol2inv(chol(mcmc[[ch]]$InvSigma[i,,]))
          sdSigmas[1+(i-B-1)/thin,] <- sqrt(diag(Sigmas[1+(i-B-1)/thin,,]))
          corSigmas[1+(i-B-1)/thin,,] <- diag(1/sdSigmas[1+(i-B-1)/thin,]) %*% Sigmas[1+(i-B-1)/thin,,] %*% diag(1/sdSigmas[1+(i-B-1)/thin,])
          detSigma[1+(i-B-1)/thin] <- det(Sigmas[1+(i-B-1)/thin,,])
        }
        for(i in 1:dim(mcmc[[ch]][["InvSigma"]])[2]){
          for(j in i:dim(mcmc[[ch]][["InvSigma"]])[2]){
            addtodatach <- as.data.frame(Sigmas[,i,j])
            colnames(addtodatach) <- paste0("Sigma[",rannames[i],",",rannames[j],"]")
            datach <- cbind(datach, addtodatach)
          }
        }
        for(i in 1:(dim(mcmc[[ch]][["InvSigma"]])[2]-1)){
          for(j in (i+1):dim(mcmc[[ch]][["InvSigma"]])[2]){
            addtodatach <- as.data.frame(corSigmas[,i,j])
            colnames(addtodatach) <- paste0("corSigma[",rannames[i],",",rannames[j],"]")
            datach <- cbind(datach, addtodatach)
          }
        }
        for(i in 1:dim(mcmc[[ch]][["InvSigma"]])[2]){
          addtodatach <- as.data.frame(sdSigmas[,i])
          colnames(addtodatach) <- paste0("sdSigma[",rannames[i],"]")
          datach <- cbind(datach, addtodatach)
        }
        addtodatach <- data.frame(detSigma = detSigma)
        datach <- cbind(datach, addtodatach)
      }
    }
    # detInvSigma, detInvQ
    for(matdet in c("detInvSigma", "detInvQ")){
      if(is.element(matdet, selected)){
        if(mcmc$spec[substring(matdet,4)]){
          for(k in 1:mcmc$K){
            dets <- numeric(length(seq(B+1, B+M, by = thin)))
            for(i in seq(B+1, B+M, by = thin)){
              dets[i] <- det(mcmc[[ch]][[substring(matdet,4)]][[k]][i,,])
            }
            addtodatach <- as.data.frame(dets)
            colnames(addtodatach) <- paste0(matdet, "(",k,")")
            datach <- cbind(datach, addtodatach)
          }
        }else{
          dets <- numeric(length(seq(B+1, B+M, by = thin)))
          for(i in seq(B+1, B+M, by = thin)){
            dets[i] <- det(mcmc[[ch]][[substring(matdet,4)]][i,,])
          }
          addtodatach <- as.data.frame(dets)
          colnames(addtodatach) <- matdet
          datach <- cbind(datach, addtodatach)
        }
      }
    }
    # random effects
    if(is.element("b",selected)){
      for(j in 1:length(unique(Y[,Id]))){
        addtodatach <- as.data.frame(mcmc[[ch]]$b[seq(B+1,B+M,by=thin),j,])
        colnames(addtodatach) <- paste0("b_[",j,",",rannames,"]")
        datach <- cbind(datach, addtodatach)
      }
    }
    #AllData <- rbind(AllData, transform(datach, chain = ch))
    datach <- datach[,-1]
    AllData[[ch]] <- as.mcmc(datach)
  } # end of for ch in ChainsToUse
  return(AllData)
}


FindAllData <- function(mcmc, B, M, what, kspec, dimspec, yspec, thin, ChainsToPlot){
  AllData <- data.frame()
  if(missing(B)){B <- mcmc$B}
  if(missing(M)){M <- mcmc$BM - B}
  
  ### All is saved in structured list
  if(mcmc$howsave == "list"){
    
    if((is.element(what, names(mcmc$spec)) & mcmc$spec[what])
       | (is.element(what, c("detInvSigma", "detSigma", "Sigma", "sdSigma", "corSigma")) & mcmc$spec["InvSigma"])
       | (is.element(what, c("detInvQ", "detQ", "Q")) & mcmc$spec["InvQ"])){
      if (is.list(mcmc[[1]][[what]][[kspec]])){
        if(missing(dimspec)){
          for(ch in ChainsToPlot){
            AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[kspec]][[yspec]][seq(B+1, B+M, by = thin)], chain = ch))
          }
        }else{
          if (length(dimspec) == 1){
            for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[kspec]][[yspec]][seq(B+1, B+M, by = thin),dimspec[1]], chain = ch))
            }
          }
          if (length(dimspec)==2){
            for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[kspec]][[yspec]][seq(B+1, B+M, by = thin),dimspec[1],dimspec[2]], chain = ch))
            }
          }
        }
        
      }else{
        if(missing(dimspec)){
          for(ch in ChainsToPlot){
            AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[kspec]][seq(B+1, B+M, by = thin)], chain = ch))
          }
        }else{
          if (length(dimspec) == 1){
            for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[kspec]][seq(B+1, B+M, by = thin),dimspec[1]], chain = ch))
            }
          }
          if (length(dimspec)==2){
            for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[kspec]][seq(B+1, B+M, by = thin),dimspec[1],dimspec[2]], chain = ch))
            }
          }
        } # end of else of if missing dimspec
        
      } # end of else of if is kspecific
    }else{
      if (is.list(mcmc[[1]][[what]])){
        if(missing(dimspec)){
          for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[yspec]][seq(B+1, B+M, by = thin)], chain = ch))
            }
        }else{
          if (length(dimspec) == 1){
            for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[yspec]][seq(B+1, B+M, by = thin),dimspec[1]], chain = ch))
            }
          }
          if (length(dimspec)==2){
            for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][[yspec]][seq(B+1, B+M, by = thin),dimspec[1],dimspec[2]], chain = ch))
            }
          }
        }
        
      }else{
        if(missing(dimspec)){
          for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][seq(B+1, B+M, by = thin)], chain = ch))
            }
        }else{
          if (length(dimspec) == 1){
            for(ch in ChainsToPlot){
                AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][seq(B+1, B+M, by = thin),dimspec[1]], chain = ch))
            }
          }
          if (length(dimspec)==2){
            for(ch in ChainsToPlot){
              AllData <- rbind(AllData, data.frame(y=mcmc[[ch]][[what]][seq(B+1, B+M, by = thin),dimspec[1],dimspec[2]], chain = ch))
            }
          }
        }
      }
    } 
  }# end of if mcmc$howsave == list
  
  ### All is saved in large matrix
  if(mcmc$howsave == "matrix"){
    #mcmc, B, M, what, kspec, dimspec, yspec, thin, ChainsToPlot
    colname <- what
    
    if(!missing(yspec)){
      colname <- paste(colname, yspec, sep = "_")
    }
    
    if(!missing(kspec)){
      colname <- paste0(colname, "(", kspec, ")")
    }
    
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        colname <- paste0(colname,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        colname <- paste0(colname,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
    
    if(!is.element(colname, colnames(mcmc$all))){
      stop(paste("Wrong specification, column ", colname, " does not exist in mcmc$all!"))
    }
    
    AllData <- data.frame()
    
    for(ch in ChainsToPlot){
      AllData <- rbind(AllData, data.frame(y=mcmc$all[seq(B+1, B+M, by = thin) + (ch-1)*mcmc$BM, colname], 
                                           chain = ch))
    }
    
  }
  return(AllData)
}


library("colorspace") 

plot.traceplots <- function(B, M, mcmc, AllData,
                            what, kspec, dimspec, yspec,
                            thin = 1,
                            ChainsToPlot,
                            COL, labcex = 1, whatLAB){
  par(mar = c(3.5,3.5,0.8,0.8))
  if (missing(ChainsToPlot)){ChainsToPlot <- c(1:min(3,mcmc$Nchains))}
  if (missing(B)){B <- mcmc$B}
  if (missing(M)){M <- mcmc$BM - B}
  #if (missing(kspec)){kspec <- 1}
  if (!missing(yspec)){if(!is.element(yspec, c(mcmc$Nums, mcmc$Ords, mcmc$Bins))){stop("Unexisting response variable specified.")}}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  #if (missing(yspec)){yspec <- ""}
  if (missing(AllData)){
    AllData <- FindAllData(mcmc = mcmc, M = M, B = B, thin = thin, 
                           what = what, kspec = kspec, dimspec = dimspec, yspec = yspec,
                           ChainsToPlot = ChainsToPlot)
  }
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(kspec)){
      whatLAB <- paste0(whatLAB, "(", kspec, ")")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  
  Max <- max(AllData$y)
  Min <- min(AllData$y)
  
  if(missing(COL)){
    if(mcmc$Nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$Nchains, c = 80, l =70)
    }
  }
  
  plot(AllData$y[1:M] ~ c((B+1):(B+M)), 
       ylim = c(Min-0.05*(Max-Min), Max+0.05*(Max-Min)),
       ylab = "", xlab = "", type= "n")
  mtext(text = "Iteration", side = 1, line = 2.5, cex = labcex)
  mtext(text = whatLAB, side = 2, line = 2.5, cex = labcex)
  for(ch in ChainsToPlot){
    lines(AllData$y[AllData$chain == ch] ~ c(seq(B+1, B+M, by = thin)),
          col = COL[ChainsToPlot[ch]])
  }
  # legend
  #legend("top", legend = ChainsToPlot, ncol = length(ChainsToPlot),
  #       col = COL[ChainsToPlot], lty =1, bty = "n", title = "Chain", cex = 0.8)
}



plot.kerneldensity <- function(B, M, mcmc, AllData,
                               what, kspec, dimspec, yspec,
                               thin = 1,
                               ChainsToPlot,
                               COL, labcex = 1, whatLAB){
  par(mar = c(3.5,3.5,0.8,0.8))
  if (missing(ChainsToPlot)){ChainsToPlot <- c(1:min(3,mcmc$Nchains))}
  if (missing(B)){B <- mcmc$B}
  if (missing(M)){M <- mcmc$BM - B}
  #if (missing(kspec)){kspec <- 1}
  if (!missing(yspec)){if(!is.element(yspec, c(mcmc$Nums, mcmc$Ords, mcmc$Bins))){stop("Unexisting response variable specified.")}}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  #if (missing(yspec)){yspec <- ""}
  if (missing(AllData)){
    AllData <- FindAllData(mcmc = mcmc, M = M, B = B, thin = thin, 
                           what = what, kspec = kspec, dimspec = dimspec, yspec = yspec,
                           ChainsToPlot = ChainsToPlot)
  }
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(kspec)){
      whatLAB <- paste0(whatLAB, "(", kspec, ")")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  if(missing(COL)){
    if(mcmc$Nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$Nchains, c = 80, l =70)
    }
  }
  
  kernest <- data.frame()
  for(ch in ChainsToPlot){
    dd <- density(AllData$y[AllData$chain == ch])
    kernest <- rbind(kernest, data.frame(x=dd$x, y=dd$y, chain = ch))
  }
  plot(c(min(kernest$y),max(kernest$y)) ~ quantile(kernest$x, probs = c(0.03,0.97)), type = "n",
       xlab = "", ylab = "")
  mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
  mtext(text = "Density", side = 2, line = 2.5, cex = labcex)
  for(ch in ChainsToPlot){
    lines(kernest$y[kernest$chain==ch] ~ kernest$x[kernest$chain == ch],
          col = COL[ch], lty = 1, lwd = 2)
  }
  #theme_set(theme_classic())
  #options(scipen = 10000)
  #require(scales)
  #g <- ggplot(AllData, aes(y)) + 
  #    geom_density(aes(fill=factor(chain)), alpha=0.7) + 
  #    labs(x = whatLAB,
  #         y = "Density",
  #         fill = "Chain") + 
  #    scale_fill_manual(values = COL[ChainsToPlot]) 
  #g
  legend("topright", legend = ChainsToPlot,
         col = COL[ChainsToPlot], lty =1, bty = "n", title = "Chain", cex = 0.8)
  
}


plot.ACF <- function(B, M, mcmc, AllData,
                     what, kspec, dimspec, yspec,
                     thin = 1, MaxLag,
                     ChainsToPlot,
                     COL, posun.width = 0.3, labcex = 1, whatLAB){
  par(mar = c(3.5,3.5,0.8,0.8))
  if (missing(ChainsToPlot)){ChainsToPlot <- c(1:min(3,mcmc$Nchains))}
  if (missing(B)){B <- mcmc$B}
  if (missing(M)){M <- mcmc$BM - B}
  #if (missing(kspec)){kspec <- 1}
  if (!missing(yspec)){if(!is.element(yspec, c(mcmc$Nums, mcmc$Ords, mcmc$Bins))){stop("Unexisting response variable specified.")}}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  #if (missing(yspec)){yspec <- ""}
  if (missing(MaxLag)){MaxLag <- min(30, M-1)}
  if (MaxLag >= M){MaxLag <- M-1}
  if (missing(AllData)){
    AllData <- FindAllData(mcmc = mcmc, M = M, B = B, thin = thin, 
                           what = what, kspec = kspec, dimspec = dimspec, yspec = yspec,
                           ChainsToPlot = ChainsToPlot)
  }
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(kspec)){
      whatLAB <- paste0(whatLAB, "(", kspec, ")")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  
  if(missing(COL)){
    if(mcmc$Nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$Nchains, c = 80, l =70)
    }
  }
  
  plot(c(-0.2,1) ~ c(0,MaxLag), 
       ylim = c(-0.2,1),
       ylab = "", xlab = "", type= "n")
  mtext(text = "Lag", side = 1, line = 2.5, cex = labcex)
  mtext(text = "ACF", side = 2, line = 2.5, cex = labcex)
  #mtext(whatLAB, side = 3, line = 0.5)
  pom <- seq(-posun.width, posun.width, length.out = length(ChainsToPlot)+1)
  posun <- (pom[1:length(ChainsToPlot)]+pom[2:(1+length(ChainsToPlot))])/2
  chcount <- 0
  for(ch in ChainsToPlot){
    chcount <- chcount+1
    pomACF <- acf(AllData$y[AllData$chain == ch], lag.max = MaxLag, plot = F)$acf[,1,1]
    lines(pomACF ~ c(0:MaxLag + posun[chcount]),
          col = COL[ch], lwd = 2)
    segments(x0 = c(0:MaxLag)+posun[chcount], 
             y0 = 0, y1 = pomACF,
             col = COL[ch], lty = 1, lwd = 1)
    abline(h = 0, col = "grey", lty = 2)
  }
}


plot.ECDF <- function(B, M, mcmc, AllData,
                      what, kspec, dimspec, yspec,
                      thin = 1, 
                      ChainsToPlot,
                      COL, labcex = 1,whatLAB){
  par(mar = c(3.5,3.5,0.8,0.8))
  if (missing(ChainsToPlot)){ChainsToPlot <- c(1:min(3,mcmc$Nchains))}
  if (missing(B)){B <- mcmc$B}
  if (missing(M)){M <- mcmc$BM - B}
  #if (missing(kspec)){kspec <- 1}
  if (!missing(yspec)){if(!is.element(yspec, c(mcmc$Nums, mcmc$Ords, mcmc$Bins))){stop("Unexisting response variable specified.")}}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  #if (missing(yspec)){yspec <- ""}
  if (missing(AllData)){
    AllData <- FindAllData(mcmc = mcmc, M = M, B = B, thin = thin, 
                           what = what, kspec = kspec, dimspec = dimspec, yspec = yspec,
                           ChainsToPlot = ChainsToPlot)
  }
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(kspec)){
      whatLAB <- paste0(whatLAB, "(", kspec, ")")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  #if(missing(whatLAB)){
  #  if(is.element(what, names(mcmc$spec)) & mcmc$spec[what]){
  #    whatLAB <- paste0(what, 
  #                      ifelse(is.list(mcmc[[1]][[what]][[kspec]]), ifelse(yspec=="",paste0("_",yspec), ""), ""),
  #                      "(",kspec,")",
  #                      ifelse(is.na(dimspec[1]),"",
  #                             paste0("[",paste(dimspec,collapse=","),"]")))
  #  }else{
  #    whatLAB <- paste0(what, 
  #                      ifelse(is.list(mcmc[[1]][[what]]), ifelse(yspec=="",paste0("_",yspec), ""), ""),
  #                      ifelse(is.na(dimspec[1]),"",
  #                             paste0("[",paste(dimspec,collapse=","),"]")))
  #  }
  #}
  
  if(missing(COL)){
    if(mcmc$Nchains == 1){
      COL <- "grey"
    }else{
      COL <- diverge_hcl(mcmc$Nchains, c = 80, l =70)
    }
  }
  #if (missing(COL)){COL <- topo.colors(mcmc$Nchains, alpha = 1)}
  
  plot(x = 1, y = 1, type = "n", xlim = c(min(AllData$y), max(AllData$y)), ylim = c(0,1),
       xlab = "", ylab = "")
  mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
  mtext(text = "ECDF", side = 2, line = 2.5, cex = labcex)
  for(ch in ChainsToPlot){
    plot(ecdf(AllData$y[AllData$chain == ch]), 
         col = COL[ch], add = T, verticals = T, pch = NA)
  }
}


plot.all <- function(B, M, mcmc,
                     what, kspec, dimspec, yspec,
                     thin = 1, MaxLag,
                     ChainsToPlot, trueval, 
                     COL, posun.width = 0.3, whatLAB){
  par(mfrow = c(2,2))
  if (missing(ChainsToPlot)){ChainsToPlot <- c(1:min(3,mcmc$Nchains))}
  if (missing(B)){B <- mcmc$B}
  if (missing(M)){M <- mcmc$BM - B}
  #if (missing(kspec)){kspec <- 1}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  if (!missing(yspec)){
    if (!is.element(yspec, c(mcmc$Nums, mcmc$Ords, mcmc$Bins))){
      stop("Unexisting response variable specified.")
    }
  }
  if (missing(MaxLag)){MaxLag <- min(30, M-1)}
  if (MaxLag >= M){MaxLag <- M-1}
  AllData <- FindAllData(mcmc = mcmc, M = M, B = B, thin = thin, 
                         what = what, kspec = kspec, dimspec = dimspec, yspec = yspec,
                         ChainsToPlot = ChainsToPlot)
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(!missing(kspec)){
      whatLAB <- paste0(whatLAB, "(", kspec, ")")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  
  plot.traceplots(B=B, M=M, mcmc=mcmc, AllData = AllData,
                  what=what, kspec=kspec, dimspec=dimspec, yspec=yspec,
                  thin=thin, ChainsToPlot=ChainsToPlot, COL=COL, labcex = 0.8)
  if(!missing(trueval)){
    abline(h = trueval, col = "seagreen", lty = 2)
  }
  plot.ECDF(B=B, M=M, mcmc=mcmc, AllData = AllData,
            what=what, kspec=kspec, dimspec=dimspec, yspec=yspec,
            thin=thin, ChainsToPlot=ChainsToPlot, COL=COL, labcex = 0.8)
  if(!missing(trueval)){
    abline(v = trueval, col = "seagreen", lty = 2)
  }
  plot.kerneldensity(B=B, M=M, mcmc=mcmc, AllData = AllData,
                     what=what, kspec=kspec, dimspec=dimspec, yspec=yspec,
                     thin=thin, ChainsToPlot=ChainsToPlot, COL=COL, labcex = 0.8)
  if(!missing(trueval)){
    abline(v = trueval, col = "seagreen", lty = 2)
  }
  plot.ACF(B=B, M=M, mcmc=mcmc, AllData = AllData,
           what=what, kspec=kspec, dimspec=dimspec, yspec=yspec,
           thin=thin, ChainsToPlot=ChainsToPlot, COL=COL, labcex = 0.8,
           MaxLag = MaxLag, posun.width = posun.width)
}


### Class-comparison - density and ECDF
plot.classes <- function(B, M, mcmc, AllData,
                         what, dimspec, yspec,
                         thin = 1, parmfrow = T,
                         ChainsToPlot, doKern = T, doECDF = T,
                         COL, labcex = 0.8, whatLAB, wherelegend = "topright"){
  if(parmfrow){par(mfrow = c(1,sum(c(doKern, doECDF))))}
  if(is.element(what,names(mcmc$spec))){
    if (!mcmc$spec[what]){stop("This parameter wasn't considered as class-specific.")}
  }
  if (missing(ChainsToPlot)){ChainsToPlot <- c(1:min(3,mcmc$Nchains))}
  if (missing(B)){B <- mcmc$B}
  if (missing(M)){M <- mcmc$BM - B}
  #if (missing(yspec)){yspec <- ifelse(is.null(mcmc$Nums), mcmc$Ords[1], mcmc$Nums[1])}
  if (!missing(yspec)){
    if (!is.element(yspec, c(mcmc$Nums, mcmc$Ords, mcmc$Bins))){
      stop("Unexisting response variable specified.")
    }
  }
  if (missing(AllData)){
    AllData <- data.frame()
    for(k in 1:mcmc$K){
      pomData <- FindAllData(mcmc = mcmc, M = M, B = B, thin = thin, 
                             what = what, kspec = k, dimspec = dimspec, yspec = yspec,
                             ChainsToPlot = ChainsToPlot)
      colnames(pomData) <- c(paste0("y(",k,")"), "chain")
      if(is.element("chain", colnames(AllData))){AllData <- cbind(AllData, pomData)
      }else{AllData <- pomData}
    }
  }
  
  if (missing(whatLAB)){
    whatLAB <- what
    
    if(!missing(yspec)){
      whatLAB <- paste(whatLAB, yspec, sep = "_")
    }
    
    if(missing(dimspec)){
      # no dimension supplied = no need for supplying -> 1-dim parameter
    }else{
      if(length(dimspec) == 1){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],"]")
      }
      if(length(dimspec) == 2){
        whatLAB <- paste0(whatLAB,"[",dimspec[1],",",dimspec[2],"]")
      }
    }
  }
  if (missing(COL)){COL <- rainbow_hcl(mcmc$K, c = 80, l =70)}
  
  # Density
  if(doKern){
    kernest <- data.frame()
    for(k in 1:mcmc$K){
      dd <- density(AllData[is.element(AllData$chain, ChainsToPlot), paste0("y(",k,")")])
      kernest <- rbind(kernest, data.frame(x=dd$x, y=dd$y, class = k))
    }
    plot(c(min(kernest$y),max(kernest$y)) ~ quantile(kernest$x, probs = c(0.03,0.97)), type = "n",
         xlab = "", ylab = "")
    mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
    mtext(text = "Density", side = 2, line = 2.5, cex = labcex)
    for(k in 1:mcmc$K){
      lines(kernest$y[kernest$class==k] ~ kernest$x[kernest$class==k],
            col = COL[k], lty = 1, lwd = 2)
    }
    legend(wherelegend, legend = 1:mcmc$K,
           col = COL, lty =1, bty = "n", title = "Class", cex = 0.8)
  }
  # ECDF
  if(doECDF){
    plot(x = 1, y = 1, type = "n", xlim = c(min(AllData[,grep("y",colnames(AllData))]), max(AllData[,grep("y",colnames(AllData))])), ylim = c(0,1),
         xlab = "", ylab = "")
    mtext(text = whatLAB, side = 1, line = 2.5, cex = labcex)
    mtext(text = "ECDF", side = 2, line = 2.5, cex = labcex)
    for(k in 1:mcmc$K){
      plot(ecdf(AllData[,paste0("y(", k,")")]), 
           col = COL[k], add = T, verticals = T, pch = NA)
    }
  }
}
