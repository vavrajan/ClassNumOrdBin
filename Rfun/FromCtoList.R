# Function to convert long vector of numbers coming from C function
# into structured list 
FromCtoList <- function(values,       # vector of values in given order
                        isspec = F,   # should [[k]] be done?
                        K = 2,        # number of classes
                        isy = F,      # should [[y]] be done?
                        ynums = 0,    # should it include numeric variables?
                        yords = 0,    # should it include ordinal variables?
                        ybins = 0,    # should it include binary variables?
                        BM = 1000,    # number of generated states 
                        d1 = 0,       # first dimension (0 means just a point state)
                        d2 = 0,       # second dimension (0 means just a vector)
                        # d1 and d2 may be dependent on y value...
                        # then d1 and d2 are vectors of possible dimensions for valid Y's
                        BYROW = T,    # by row or by column?
                        sym = 0,      # if D = 2, is it symmetric matrix?
                        diag = 0,     # If symmetrical, is diagonal stored?
                        diagval = 1,  # If sym and not diag, what is the diag value?
                        D = 0)        # dimension of converted parameter
                                      # D = 0 --> parameter is just one number
                                      #           (can be stored as a vector)
                                      # D = 1 --> parameter is a vector
                                      #           (stored as matrix, values in columns)
                                      # D = 2 --> parameter is a matrix 
                                      #           (stored as array, stored as [1,.,.])
{
  Ys <- c()
  if(ynums){Ys <- c(Ys, Nums)}
  if(yords){Ys <- c(Ys, Ords)}
  if(ybins){Ys <- c(Ys, Bins)}
  
  if(isspec){
    RET <- list()
    if(isy){
      
      if(D==0){
        for(k in 1:K){
          ## option with [[y]]
          #RET[[k]] <- list()
          #for(y in 1:length(Ys)){
            #RET[[k]][[Ys[y]]] <- values[seq(from = (k-1)*length(Ys)+y, 
            #                                to = K*length(Ys)*BM, 
            #                                by = K*length(Ys))]
            # }
          
          ## better --> no need for [[y]] --> can be stored in matrix with columns names Ys
          RET[[k]] <- matrix(0, nrow = BM, ncol = length(Ys))
          colnames(RET[[k]]) <- Ys
          for(y in 1:length(Ys)){
            RET[[k]][, Ys[y]] = values[seq(from = (k-1)*length(Ys)+y, 
                                           to = K*length(Ys)*BM, 
                                           by = K*length(Ys))]
          }
            
        }
      } # end of D==0
      
      if(D==1){
        for(k in 1:K){
          RET[[k]] <- list()
          dimy = 0
          for(y in 1:length(Ys)){
            RET[[k]][[Ys[y]]] <- matrix(0, nrow = BM, ncol = d1[Ys[y]])
            for(i in 1:BM){
              RET[[k]][[Ys[y]]][i, ] <- values[(i-1)*K*sum(d1) + (k-1)*sum(d1) + dimy + 1:d1[Ys[y]]]
            }
            dimy = dimy + d1[Ys[y]]
          }
        }
      } # end of D==1
      
      if(D==2){
        for(k in 1:K){
          RET[[k]] <- list()
          if(sym){
            if(diag){
              dimtot = sum(d1*(d1+1)/2)
            }else{
              dimtot = sum((d1-1)*d1/2)
            }
          }else{
            dimtot = sum(d1*d2)
          }
          dimy = 0
          for(y in 1:length(Ys)){
            RET[[k]][[Ys[y]]] <- array(0, dim = c(BM, d1[Ys[y]], d2[Ys[y]]))
            for(i in 1:BM){
              if(sym){
                # matrix is symmetrical and only upper-right triangle is stored
                # d1 = d2
                if(diag){
                  dd = d1[Ys[y]]*(d1[Ys[y]]+1)/2 # number of elements 
                  for(row in 1:d1[Ys[y]]){
                    for(col in row:d1[Ys[y]]){
                      RET[[k]][[Ys[y]]][i, row, col] <- RET[[k]][[Ys[y]]][i, col, row] <-
                        values[(i-1)*K*dimtot + (k-1)*dimtot + dimy + row + (col-1)*col/2]
                    }
                  }
                }else{
                  dd = (d1[Ys[y]]-1)*d1[Ys[y]]/2 # number of elements
                  # non-diagonal elements
                  for(row in 1:(d1[Ys[y]]-1)){
                    for(col in (row+1):d1[Ys[y]]){
                      RET[[k]][[Ys[y]]][i, row, col] <- RET[[k]][[Ys[y]]][i, col, row] <-
                        values[(i-1)*K*dimtot + (k-1)*dimtot + dimy + row + (col-2)*(col-1)/2]
                    }
                  }
                  # diagonal elements
                  for(row in 1:d1[Ys[y]]){
                    RET[[k]][[Ys[y]]][i, row, row] <- diagval
                  }
                } # end of else of diag
              }else{
                # matrix is general rectangular --> all is stored
                dd = d1[Ys[y]]*d2[Ys[y]]
                RET[[k]][[Ys[y]]][i,,] <- matrix(values[(i-1)*K*dimtot + (k-1)*dimtot + dimy + 1:dd],
                                                 nrow = d1[Ys[y]],
                                                 ncol = d2[Ys[y]],
                                                 byrow = BYROW)
              } # end of else of sym
            }
            dimy = dimy + dd
          }
        }
      } # end of D==2
      
    }else{
      # is class-specific, however not y-specific
      
      if(D==0){
        RET <- matrix(0, nrow = BM, ncol = K)
        for(k in 1:K){
          # no need to do [[k]] --> can be stored in matrix
          # column corresponds to k
          RET[, k] = values[seq(from = k,
                                to = K*BM,
                                by = K)]
          
        }
      } # end of D==0
      
      if(D==1){
        for(k in 1:K){
          RET[[k]] <- matrix(0, nrow = BM, ncol = d1)
          for(i in 1:BM){
            RET[[k]][i, ] <- values[(i-1)*K*d1 + (k-1)*d1 + 1:d1]          
          }
        }
      } # end of D==1
      
      if(D==2){
        for(k in 1:K){
          RET[[k]] <- array(0, dim = c(BM, d1, d2))
          for(i in 1:BM){
            if(sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(diag){
                dd = d1*(d1+1)/2 # number of elements 
                for(row in 1:d1){
                  for(col in row:d1){
                    RET[[k]][i, row, col] <- RET[[k]][i, col, row] <-
                      values[(i-1)*K*dd + (k-1)*dd + row + (col-1)*col/2]
                  }
                }
              }else{
                dd = (d1-1)*d1/2 # number of elements
                # non-diagonal elements
                for(row in 1:(d1-1)){
                  for(col in (row+1):d1){
                    RET[[k]][i, row, col] <- RET[[k]][i, col, row] <-
                      values[(i-1)*K*dd + (k-1)*dd + row + (col-2)*(col-1)/2]
                  }
                }
                # diagonal elements
                for(row in 1:d1){
                  RET[[k]][i, row, row] <- diagval
                }
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              RET[[k]][i,,] <- matrix(values[(i-1)*K*d1*d2 + (k-1)*d1*d2 + 1:(d1*d2)],
                                               nrow = d1,
                                               ncol = d2,
                                               byrow = BYROW)
            } # end of else of sym
          }
        }
      } # end of D==2
      
    } # end of else of isy
  }else{
    # not class-specific
    if(isy){
      # not class-specific, but still y-specific
      if(D==0){
          ## option with [[y]]
          #RET <- list()
          #for(y in 1:length(Ys)){
          #RET[[Ys[y]]] <- values[seq(from = y, 
          #                                to = length(Ys)*BM, 
          #                                by = length(Ys))]
          # }
          
          ## better --> no need for [[y]] --> can be stored in matrix with columns names Ys
          RET <- matrix(0, nrow = BM, ncol = length(Ys))
          colnames(RET) <- Ys
          for(y in 1:length(Ys)){
            RET[, Ys[y]] = values[seq(from = y,
                                      to = length(Ys)*BM,
                                      by = length(Ys))]
          }
          
      } # end of D==0
      
      if(D==1){
        RET <- list()
        dimtot = sum(d1)
        dimy = 0
        for(y in 1:length(Ys)){
          RET[[Ys[y]]] <- matrix(0, nrow = BM, ncol = d1[Ys[y]])
          for(i in 1:BM){
            RET[[Ys[y]]][i, ] <- values[(i-1)*dimtot + dimy + 1:d1[Ys[y]]]
          }
          dimy = dimy + d1[Ys[y]]
        }
      } # end of D==1
      
      if(D==2){
        RET <- list()
        if(sym){
          if(diag){
            dimtot = sum(d1*(d1+1)/2)
          }else{
            dimtot = sum((d1-1)*d1/2)
          }
        }else{
          dimtot = sum(d1*d2)
        }
        dimy = 0
        for(y in 1:length(Ys)){
          RET[[Ys[y]]] <- array(0, dim = c(BM, d1[Ys[y]], d2[Ys[y]]))
          for(i in 1:BM){
            if(sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(diag){
                dd = d1[Ys[y]]*(d1[Ys[y]]+1)/2 # number of elements 
                for(row in 1:d1[Ys[y]]){
                  for(col in row:d1[Ys[y]]){
                    RET[[Ys[y]]][i, row, col] <- RET[[Ys[y]]][i, col, row] <-
                      values[(i-1)*dimtot + dimy + row + (col-1)*col/2]
                  }
                }
              }else{
                dd = (d1[Ys[y]]-1)*d1[Ys[y]]/2 # number of elements
                # non-diagonal elements
                for(row in 1:(d1[Ys[y]]-1)){
                  for(col in (row+1):d1[Ys[y]]){
                    RET[[Ys[y]]][i, row, col] <- RET[[Ys[y]]][i, col, row] <-
                      values[(i-1)*dimtot  + dimy + row + (col-2)*(col-1)/2]
                  }
                }
                # diagonal elements
                for(row in 1:d1[Ys[y]]){
                  RET[[Ys[y]]][i, row, row] <- diagval
                }
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              dd = d1[Ys[y]] * d2[Ys[y]]
              RET[[Ys[y]]][i,,] <- matrix(values[(i-1)*dimtot + dimy + 1:dd],
                                          nrow = d1[Ys[y]],
                                          ncol = d2[Ys[y]],
                                          byrow = BYROW)
            } # end of else of sym
          }
          dimy = dimy + dd
        }
      } # end of D==2
      
    }else{
      # is NEITHER class-specific, nor y-specific
      
      if(D==0){
        RET <- values
      } # end of D==0
      
      if(D==1){
        RET <-  matrix(values, nrow = BM, ncol = d1, byrow = BYROW)
      } # end of D==1
      
      if(D==2){
        RET <- array(0, dim = c(BM, d1, d2))
        for(i in 1:BM){
          if(sym){
            # matrix is symmetrical and only upper-right triangle is stored
            # d1 = d2
            if(diag){
              dd = d1*(d1+1)/2 # number of elements 
              for(row in 1:d1){
                for(col in row:d1){
                  RET[i, row, col] <- RET[i, col, row] <-
                    values[(i-1)*dd + row + (col-1)*col/2]
                }
              }
            }else{
              dd = (d1-1)*d1/2 # number of elements
              # non-diagonal elements
              for(row in 1:(d1-1)){
                for(col in (row+1):d1){
                  RET[i, row, col] <- RET[i, col, row] <-
                    values[(i-1)*dd + row + (col-2)*(col-1)/2]
                }
              }
              # diagonal elements
              for(row in 1:d1){
                RET[i, row, row] <- diagval
              }
            } # end of else of diag
          }else{
            # matrix is general rectangular --> all is stored
            RET[i,,] <- matrix(values[(i-1)*d1*d2 + 1:(d1*d2)],
                                    nrow = d1,
                                    ncol = d2,
                                    byrow = BYROW)
          } # end of else of sym
        }
      } # end of D==2
    } # end of else of isy
  } # end of else of isspec
  
  return(RET)
}





FromCtoList_settings <- function(values, 
                                 p,
                                 settings,
                                 yspecd1, 
                                 yspecd2)        # parameter to be listed
{
  # data frame settings consists of needed parameters:
  # save, isspec, K, isy, ynums, yords, ybins, BM, 
  # d1, d2, BYROW, sym, diag, diagval, D
  v <- c(settings[p, ])
  Ys <- c()
  if(v$ynums){Ys <- c(Ys, Nums)}
  if(v$yords){Ys <- c(Ys, Ords)}
  if(v$ybins){Ys <- c(Ys, Bins)}
  
  sumd1 <- ifelse(missing(yspecd1), v$d1*length(Ys), sum(yspecd1))
  #sumd1 <- v$d1*length(Ys)
  
  if(v$isspec){
    RET <- list()
    if(v$isy){
      
      if(v$D==0){
        for(k in 1:v$K){
          ## option with [[y]]
          #RET[[k]] <- list()
          #for(y in 1:length(Ys)){
          #RET[[k]][[Ys[y]]] <- values[seq(from = (k-1)*length(Ys)+y, 
          #                                to = K*length(Ys)*BM, 
          #                                by = K*length(Ys))]
          # }
          
          ## better --> no need for [[y]] --> can be stored in matrix with columns names Ys
          RET[[k]] <- matrix(0, nrow = v$BM, ncol = length(Ys))
          colnames(RET[[k]]) <- Ys
          for(y in 1:length(Ys)){
            RET[[k]][, Ys[y]] = values[seq(from = (k-1)*length(Ys)+y, 
                                           to = v$K*length(Ys)*v$BM, 
                                           by = v$K*length(Ys))]
          }
          
        }
      } # end of D==0
      
      if(v$D==1){
        for(k in 1:v$K){
          RET[[k]] <- list()
          dimy = 0
          for(y in 1:length(Ys)){
            d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
            #d1 = v$d1
            RET[[k]][[Ys[y]]] <- matrix(0, nrow = v$BM, ncol = d1)
            for(i in 1:v$BM){
              RET[[k]][[Ys[y]]][i, ] <- values[(i-1)*v$K*sumd1 + (k-1)*sumd1 + dimy + 1:d1]
            }
            dimy = dimy + d1
          }
        }
      } # end of D==1
      
      if(v$D==2){
        for(k in 1:v$K){
          RET[[k]] <- list()
          if(v$sym){
            if(v$diag){
              dimtot = ifelse(missing(yspecd1), length(Ys)*v$d1*(v$d1+1)/2,sum(yspecd1*(yspecd1+1)/2))
            }else{
              dimtot = ifelse(missing(yspecd1), length(Ys)*(v$d1-1)*v$d1/2,sum((yspecd1-1)*yspecd1/2))
            }
          }else{
            dimtot = ifelse(missing(yspecd1)|missing(yspecd2), length(Ys)*v$d1*v$d2, sum(yspecd1*yspecd2))
          }
          dimy = 0
          for(y in 1:length(Ys)){
            d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
            d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
            RET[[k]][[Ys[y]]] <- array(0, dim = c(v$BM, d1, d2))
            for(i in 1:v$BM){
              if(v$sym){
                # matrix is symmetrical and only upper-right triangle is stored
                # d1 = d2
                if(v$diag){
                  dd = d1*(d1+1)/2 # number of elements 
                  for(row in 1:d1){
                    for(col in row:d1){
                      RET[[k]][[Ys[y]]][i, row, col] <- RET[[k]][[Ys[y]]][i, col, row] <-
                        values[(i-1)*v$K*dimtot + (k-1)*dimtot + dimy + row + (col-1)*col/2]
                    }
                  }
                }else{
                  dd = (d1-1)*d1/2 # number of elements
                  # non-diagonal elements
                  for(row in 1:(d1-1)){
                    for(col in (row+1):d1){
                      RET[[k]][[Ys[y]]][i, row, col] <- RET[[k]][[Ys[y]]][i, col, row] <-
                        values[(i-1)*v$K*dimtot + (k-1)*dimtot + dimy + row + (col-2)*(col-1)/2]
                    }
                  }
                  # diagonal elements
                  for(row in 1:d1){
                    RET[[k]][[Ys[y]]][i, row, row] <- v$diagval
                  }
                } # end of else of diag
              }else{
                # matrix is general rectangular --> all is stored
                dd = d1*d2
                RET[[k]][[Ys[y]]][i,,] <- matrix(values[(i-1)*v$K*dimtot + (k-1)*dimtot + dimy + 1:dd],
                                                 nrow = d1,
                                                 ncol = d2,
                                                 byrow = v$BYROW)
              } # end of else of sym
            }
            dimy = dimy + dd
          }
        }
      } # end of D==2
      
    }else{
      # is class-specific, however not y-specific
      
      if(v$D==0){
        RET <- matrix(0, nrow = v$BM, ncol = v$K)
        for(k in 1:v$K){
          # no need to do [[k]] --> can be stored in matrix
          # column corresponds to k
          RET[, k] = values[seq(from = k,
                                to = v$K*v$BM,
                                by = v$K)]
          
        }
      } # end of D==0
      
      if(v$D==1){
        for(k in 1:v$K){
          RET[[k]] <- matrix(0, nrow = v$BM, ncol = v$d1)
          for(i in 1:v$BM){
            RET[[k]][i, ] <- values[(i-1)*v$K*v$d1 + (k-1)*v$d1 + 1:v$d1]          
          }
        }
      } # end of D==1
      
      if(v$D==2){
        for(k in 1:v$K){
          RET[[k]] <- array(0, dim = c(v$BM, v$d1, v$d2))
          for(i in 1:v$BM){
            if(v$sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(v$diag){
                dd = v$d1*(v$d1+1)/2 # number of elements 
                for(row in 1:v$d1){
                  for(col in row:v$d1){
                    RET[[k]][i, row, col] <- RET[[k]][i, col, row] <-
                      values[(i-1)*v$K*dd + (k-1)*dd + row + (col-1)*col/2]
                  }
                }
              }else{
                dd = (v$d1-1)*v$d1/2 # number of elements
                # non-diagonal elements
                for(row in 1:(v$d1-1)){
                  for(col in (row+1):v$d1){
                    RET[[k]][i, row, col] <- RET[[k]][i, col, row] <-
                      values[(i-1)*v$K*dd + (k-1)*dd + row + (col-2)*(col-1)/2]
                  }
                }
                # diagonal elements
                for(row in 1:v$d1){
                  RET[[k]][i, row, row] <- v$diagval
                }
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              RET[[k]][i,,] <- matrix(values[(i-1)*v$K*v$d1*v$d2 + (k-1)*v$d1*v$d2 + 1:(v$d1*v$d2)],
                                      nrow = v$d1,
                                      ncol = v$d2,
                                      byrow = v$BYROW)
            } # end of else of sym
          }
        }
      } # end of D==2
      
    } # end of else of isy
  }else{
    # not class-specific
    if(v$isy){
      # not class-specific, but still y-specific
      if(v$D==0){
        ## option with [[y]]
        #RET <- list()
        #for(y in 1:length(Ys)){
        #RET[[Ys[y]]] <- values[seq(from = y, 
        #                                to = length(Ys)*BM, 
        #                                by = length(Ys))]
        # }
        
        ## better --> no need for [[y]] --> can be stored in matrix with columns names Ys
        RET <- matrix(0, nrow = v$BM, ncol = length(Ys))
        colnames(RET) <- Ys
        for(y in 1:length(Ys)){
          RET[, Ys[y]] = values[seq(from = y,
                                    to = length(Ys)*v$BM,
                                    by = length(Ys))]
        }
        
      } # end of D==0
      
      if(v$D==1){
        RET <- list()
        dimtot = ifelse(missing(yspecd1), length(Ys)*v$d1, sum(yspecd1))
        dimy = 0
        for(y in 1:length(Ys)){
          d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
          RET[[Ys[y]]] <- matrix(0, nrow = v$BM, ncol = d1)
          for(i in 1:v$BM){
            RET[[Ys[y]]][i, ] <- values[(i-1)*dimtot + dimy + 1:d1]
          }
          dimy = dimy + d1
        }
      } # end of D==1
      
      if(v$D==2){
        RET <- list()
        if(v$sym){
          if(v$diag){
            dimtot = ifelse(missing(yspecd1), length(Ys)*v$d1*(v$d1+1)/2, sum(yspecd1*(yspecd1+1)/2))
          }else{
            dimtot = ifelse(missing(yspecd1), length(Ys)*(v$d1-1)*v$d1/2, sum((yspecd1-1)*yspecd1/2))
          }
        }else{
          dimtot = ifelse(missing(yspecd1)|missing(yspecd2), length(Ys)*v$d1*v$d2, sum(yspecd1*yspecd2))
        }
        dimy = 0
        for(y in 1:length(Ys)){
          d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
          d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
          RET[[Ys[y]]] <- array(0, dim = c(v$BM, d1, d2))
          for(i in 1:v$BM){
            if(v$sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(v$diag){
                dd = d1*(d1+1)/2 # number of elements 
                for(row in 1:d1){
                  for(col in row:d1){
                    RET[[Ys[y]]][i, row, col] <- RET[[Ys[y]]][i, col, row] <-
                      values[(i-1)*dimtot + dimy + row + (col-1)*col/2]
                  }
                }
              }else{
                dd = (d1-1)*d1/2 # number of elements
                # non-diagonal elements
                for(row in 1:(d1-1)){
                  for(col in (row+1):d1){
                    RET[[Ys[y]]][i, row, col] <- RET[[Ys[y]]][i, col, row] <-
                      values[(i-1)*dimtot  + dimy + row + (col-2)*(col-1)/2]
                  }
                }
                # diagonal elements
                for(row in 1:d1){
                  RET[[Ys[y]]][i, row, row] <- v$diagval
                }
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              dd = d1 * d2
              RET[[Ys[y]]][i,,] <- matrix(values[(i-1)*dimtot + dimy + 1:dd],
                                          nrow = d1,
                                          ncol = d2,
                                          byrow = v$BYROW)
            } # end of else of sym
          }
          dimy = dimy + dd
        }
      } # end of D==2
      
    }else{
      # is NEITHER class-specific, nor y-specific
      
      if(v$D==0){
        RET <- values
      } # end of D==0
      
      if(v$D==1){
        RET <-  matrix(values, nrow = v$BM, ncol = v$d1, byrow = v$BYROW)
      } # end of D==1
      
      if(v$D==2){
        RET <- array(0, dim = c(v$BM, v$d1, v$d2))
        for(i in 1:v$BM){
          if(v$sym){
            # matrix is symmetrical and only upper-right triangle is stored
            # d1 = d2
            if(v$diag){
              dd = v$d1*(v$d1+1)/2 # number of elements 
              for(row in 1:v$d1){
                for(col in row:v$d1){
                  RET[i, row, col] <- RET[i, col, row] <-
                    values[(i-1)*dd + row + (col-1)*col/2]
                }
              }
            }else{
              dd = (v$d1-1)*v$d1/2 # number of elements
              # non-diagonal elements
              for(row in 1:(v$d1-1)){
                for(col in (row+1):v$d1){
                  RET[i, row, col] <- RET[i, col, row] <-
                    values[(i-1)*dd + row + (col-2)*(col-1)/2]
                }
              }
              # diagonal elements
              for(row in 1:v$d1){
                RET[i, row, row] <- v$diagval
              }
            } # end of else of diag
          }else{
            # matrix is general rectangular --> all is stored
            RET[i,,] <- matrix(values[(i-1)*v$d1*v$d2 + 1:(v$d1*v$d2)],
                               nrow = v$d1,
                               ncol = v$d2,
                               byrow = v$BYROW)
          } # end of else of sym
        }
      } # end of D==2
    } # end of else of isy
  } # end of else of isspec
  
  return(RET)
}
