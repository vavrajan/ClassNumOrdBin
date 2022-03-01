PermuteClusterLabels <- function(mcmc, perm){
  # perm contains list of permutations - each for one chain
  # perm[[ch]][1] contains which old cluster should be the new first one (for chain ch)
  # ... and so on
  
  if(mcmc$howsave == "matrix"){
    ## data are save in matrix mcmc$all
    for(ch in 1:mcmc$Nchains){
      # w
      if(mcmc$whatsave["w"]){
        auxw <- matrix(-1, nrow = sum(mcmc$all$chain==ch), ncol = 0)
        for(k in 1:mcmc$K){
          auxw <- cbind(auxw, mcmc$all[mcmc$all$chain==ch, paste0("w[",perm[[ch]][k],"]")])
        }
        mcmc$all[mcmc$all$chain==ch, paste0("w[",1:mcmc$K,"]")] <- auxw
      } 
      
      # pUik
      if(mcmc$whatsave["pUik"]){
        for(i in 1:mcmc$settings["U","dims"]){
          auxpUik <- matrix(-1, nrow =sum(mcmc$all$chain==ch), ncol = 0)
          auxpUi <- mcmc$all[mcmc$all$chain==ch, paste0("pUik[",i,",",1:mcmc$K, "]")]
          for(k in 1:mcmc$K){
            auxpUik <- cbind(auxpUik, auxpUi[,perm[[ch]][k]])
          }
          mcmc$all[mcmc$all$chain==ch, paste0("pUik[",i,",",1:mcmc$K, "]")] <- auxpUik
        }
      }
      
      # permute U
      if(mcmc$whatsave["U"]){
        for(i in 1:mcmc$settings["U","dims"]){
          auxU <- mcmc$all[mcmc$all$chain==ch, paste0("U[",i,"]")]+1
          auxU <- sapply(auxU, function(k){which(perm[[ch]] == k)})
          mcmc$all[mcmc$all$chain==ch, paste0("U[",i,"]")] <- auxU-1
        }
      }
      
      # cluster-specific parameters
      aux <- mcmc$all[mcmc$all$chain == ch, ]
      cnames <- colnames(aux)[grep("\\(1\\)", colnames(aux))]
      for(k in 1:mcmc$K){
        pnames <- gsub("\\(1\\)", paste0("\\(",perm[[ch]][k],"\\)"), cnames)
        knames <- gsub("\\(1\\)", paste0("\\(",k,"\\)"), cnames)
        aux[aux$chain==ch, knames] <- mcmc$all[aux$chain==ch, pnames]
      } 
      mcmc$all[mcmc$all$chain==ch, ] <- aux
    }
  }
  
  if(mcmc$howsave == "list"){
    # data are saved in lists
    for(ch in 1:mcmc$Nchains){
      # w
      if(mcmc$whatsave["w"]){
        auxw <- matrix(-1, nrow = mcmc$MB, ncol = 0)
        for(k in 1:mcmc$K){
          auxw <- cbind(auxw, mcmc[[ch]]$w[, perm[[ch]][k]])
        }
        mcmc[[ch]]$w <- auxw
      } 
      
      # pUik
      if(mcmc$whatsave["pUik"]){
        for(i in 1:mcmc$settings["U","dims"]){
          auxpUik <- matrix(-1, nrow = mcmc$MB, ncol = 0)
          auxpUi <- mcmc[[ch]]$pUik[,i,1:mcmc$K]
          for(k in 1:mcmc$K){
            auxpUik <- cbind(auxpUik, auxpUi[,perm[[ch]][k]])
          }
          mcmc[[ch]]$pUik[,i,1:mcmc$K] <- auxpUik
        }
      }
      
      # permute U
      if(mcmc$whatsave["U"]){
        for(i in 1:mcmc$settings["U","dims"]){
          auxU <- mcmc[[ch]]$U[,i]+1
          auxU <- sapply(auxU, function(k){which(perm[[ch]] == k)})
          mcmc[[ch]]$U[,i] <- auxU-1
        }
      }
      
      # cluster-specific parameters
      for(p in rownames(mcmc$settings)[as.logical(mcmc$settings$isspec)]){
        aux <- mcmc[[ch]][[p]]
        for(k in 1:mcmc$K){
          aux[[k]] <- mcmc[[ch]][[p]][[perm[[ch]][k]]]
        } 
        mcmc[[ch]][[p]] <- aux
      }
    }
  }
  
  return(mcmc)
}
