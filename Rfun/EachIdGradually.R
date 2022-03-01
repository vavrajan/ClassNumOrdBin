# Expands dataset by repetition of observations
# Each subject will have first j=1,...,n_i observations stored
# under a new Id that will be id*100+j
# observations are ordered according to variable Order
# increasingly

EachIdGradually <- function(data,
                            Id = "subject",
                            Order = "time",
                            digits_added = 2){
  # first find all individuals
  Ids <- unique(data[,Id])
  # For each individual create gradual dataset
  RES <- data.frame()
  for(id in Ids){
    subdata <- data[data[,Id] == id, ]
    subdata <- subdata[order(subdata[,Order]), ]
    nid <- dim(subdata)[1]
    for(j in 1:nid){
      newdata <- subdata[1:j, ]
      newid <- as.numeric(id)*10^digits_added + j
      newdata$newid <- rep(newid, j)
      RES <- rbind(RES, newdata)
    }
  }
  return(RES)
}
