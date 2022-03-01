GenerateData <- function(n = n, n_i = 4, 
                         K = 2, difference = "intercept", randompart = "intercept",
                         probs = probs, Betamu = Betamu, Sigma = Sigma, Eb = Eb,
                         tauNum = tauNum, gam = gam, gamBin = gamBin,
                         nY = nY, Ords = Ords, Bins = Bins){
  data <- data.frame()
  rdim <- ifelse(randompart == "both", 2*sum(nY), sum(nY))
  for (i in 1:n){
    datai <- data.frame()
    datai[1:n_i, "subject"] <- i
    newU <- sump <- 0
    u1 <- runif(1)
    while(sump < u1){
      newU <- newU+1
      sump <- sump + probs[[K]][newU]
    }
    datai[, "U"] <- newU
    datai[, "j"] <- 1:n_i
    datai[, "time"] <- sort(runif(n_i))
    b <- t(chol(Sigma[[K]][[newU]][[randompart]])) %*% rnorm(rdim,0,1) + Eb[[K]][[newU]][[randompart]][[difference]]
    
    for(y in 1:sum(nY)){
      yy <- paste0("Y",y)
      datai[,yy] <- switch(randompart,
                           intercept = b[y] + Betamu[[K]][[newU]][[randompart]][[difference]][[yy]][4] * datai[, "time"],
                           slope = Betamu[[K]][[newU]][[randompart]][[difference]][[yy]][3] + b[y] * datai[, "time"],
                           both = b[y] + b[y+sum(nY)] * datai[, "time"])
      if(is.element(randompart, c("intercept", "both"))){
        datai[, paste0("b_",yy,"_","intercept")] <- b[y] 
      }
      if(randompart == "slope"){
        datai[, paste0("b_",yy,"_","slope")] <- b[y] 
      }
      if(randompart == "both"){
        datai[, paste0("b_",yy,"_","slope")] <- b[y+sum(nY)] 
      }
    }
    
    # regressor X1 ~ Alt(0.5) - if time-invariant for each subject
    datai[, "X1"] <- rbinom(1, 1, 0.5)
    
    data <- rbind(data, datai)
  }
  
  ## regressors
  # regressor X1 ~ Alt(0.5) - if time-varying 
  #X1 <- rbinom(n*n_i, 1, prob = 0.5)
  
  # regressor X2 ~ Unif(0,1)
  data[,"X2"] <- runif(n*n_i)
  
  ## Responses
  for(y in 1:sum(nY)){
    yy <- paste0("Y",y)
    data[, yy] <- data[, yy] + rnorm(n*n_i, 0, sd = ifelse(yy==1,sqrt(1/tauNum),1)) 
    for(xx in 1:2){
      data[, yy] <- data[, yy] + apply(data, 1, function(row){
        return(as.numeric(Betamu[[K]][[as.numeric(row[["U"]])]][[randompart]][[difference]][[yy]][xx]) * as.numeric(row[paste0("X",xx)]))
      })
    }
  }
  
  for(fy in Ords){
    yy <- gsub("f","",fy)
    data[,fy] <- cut(data[,yy], breaks = c(-Inf, -1, gam, Inf), labels = c(0:2), ordered_result = T)
  }
  
  for(yy in Bins){
    fy <- paste0("f",yy)
    data[,fy] <- cut(data[,yy], breaks = c(-Inf, gamBin, Inf), labels = c(0,1), ordered_result = T)
  }
  
  return(data)
}
