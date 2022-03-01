FromMatrixtoC_settings <- function(all,
                                   chain,
                                   iterations,
                                   p,
                                   settings,
                                   yspecd1,
                                   yspecd2)
{
  suball <- all[all$chain == chain, ]
  if(missing(iterations)){iterations <- 1:dim(suball)[1]}
  
  ### We just need to take those columns that belong to parameter p
  # Every colname begins with p (--> "^")
  # be careful about "Sigma" being contained in "InvSigma",...
  # "pUik" vs. "pUik_nonb"
  # also careful about p = "b" and "beta"
  
  cols <- grep(paste0("^",p), colnames(all), value = T)
  # These cols begin with parameter p,
  # but we need to exclude those which might be onfused with other parameter
  otherp <- setdiff(rownames(settings), p)
  for(op in otherp){
    if(length(grep(p, op))>0){ # p is contained in op
      # --> all op must be deleted from this list of columns
      cols <- setdiff(cols, grep(paste0("^",op), cols, value = T))
    }
  }
  
  # matrix needs to be stored by rows --> needs to be transposed
  RET <- c(t(suball[iterations, cols]))
  return(RET)
}
