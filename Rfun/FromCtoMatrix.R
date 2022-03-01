FromCtoMatrix_settings <- function(values,
                                   p,
                                   settings,
                                   yspecd1,
                                   yspecd2,
                                   Nums, Ords, Bins)        # parameter to be listed
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
  
  RET <- as.data.frame(matrix(values, nrow = v$BM, ncol = v$dimswithK, byrow = T))
  
  # what remains is to set colnames properly
  
  
  if(v$isspec){
    if(v$isy){
      
      if(v$D==0){
        for(k in 1:v$K){
          colnames(RET)[(k-1)*length(Ys) + 1:length(Ys)] <- paste0(p,"_",Ys,"(",k,")")
        }
      } # end of D==0
      
      if(v$D==1){
        for(k in 1:v$K){
          dimy = 0
          for(y in 1:length(Ys)){
            d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
            colnames(RET)[(k-1)*sumd1 + dimy + 1:d1] <- paste0(p,"_",Ys[y],"(",k,")[",1:d1,"]")
            dimy = dimy + d1
          }
        }
      } # end of D==1
      
      if(v$D==2){
          
        for(k in 1:v$K){
          dimy = 0
          for(y in 1:length(Ys)){
            d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
            d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
            cnames <- c()
            if(v$sym){
              # matrix is symmetrical and only upper-right triangle is stored
              # d1 = d2
              if(v$diag){
                cnames <- c()
                for(col in 1:d1){
                  for(row in 1:col){
                    cnames <- c(cnames, paste0(p,"_",Ys[y],"(",k,")[",row,",",col,"]"))
                  }
                }
                dd = d1*(d1+1)/2 # number of elements
                colnames(RET)[(k-1)*v$dims + dimy + 1:dd] <- cnames
              }else{
                # only non-diagonal elements
                for(col in 2:d1){
                  for(row in 1:(col-1)){
                    cnames <- c(cnames, paste0(p,"_",Ys[y],"(",k,")[",row,",",col,"]"))
                  }
                }
                dd = (d1-1)*d1/2 # number of elements
                colnames(RET)[(k-1)*v$dims + dimy + 1:dd] <- cnames
              } # end of else of diag
            }else{
              # matrix is general rectangular --> all is stored
              if(v$BYROW){
                for(row in 1:d1){
                  cnames <- c(cnames, paste0(p,"_",Ys[y],"(",k,")[",row,",",1:d2,"]"))
                }
              }else{
                for(col in 1:d2){
                  cnames <- c(cnames, paste0(p,"_",Ys[y],"(",k,")[",1:d1,",",col,"]"))
                }
              }
              dd = d1*d2
              colnames(RET)[(k-1)*v$dims + dimy + 1:(d1*d2)] <- cnames
            } # end of else of sym
            dimy = dimy + dd
          } # end of for y
        } # end of for k
      } # end of D==2
      
    }else{
      # is class-specific, however not y-specific
      
      if(v$D==0){
        colnames(RET) <- paste0(p,"(",1:v$K,")")
      } # end of D==0
      
      if(v$D==1){
        for(k in 1:v$K){
          colnames(RET)[(k-1)*v$d1 + 1:v$d1] <- paste0(p,"(",k,")[",1:v$d1,"]")
        }
      } # end of D==1
      
      if(v$D==2){
        for(k in 1:v$K){
          cnames <- c()
          if(v$sym){
            # matrix is symmetrical and only upper-right triangle is stored
            # d1 = d2
            if(v$diag){
              for(col in 1:v$d1){
                for(row in 1:col){
                  cnames <- c(cnames, paste0(p,"(",k,")[",row,",",col,"]"))
                }
              }
              dd = v$d1*(v$d1+1)/2 # number of elements
              colnames(RET)[(k-1)*v$dims + 1:dd] <- cnames
            }else{
              # only non-diagonal elements
              for(col in 2:v$d1){
                for(row in 1:(col-1)){
                  cnames <- c(cnames, paste0(p,"(",k,")[",row,",",col,"]"))
                }
              }
              dd = (v$d1-1)*v$d1/2 # number of elements
              colnames(RET)[(k-1)*v$dims + 1:dd] <- cnames
            } # end of else of diag
          }else{
            # matrix is general rectangular --> all is stored
            if(v$BYROW){
              for(row in 1:v$d1){
                cnames <- c(cnames, paste0(p,"(",k,")[",row,",",1:v$d2,"]"))
              }
            }else{
              for(col in 1:v$d2){
                cnames <- c(cnames, paste0(p,"(",k,")[",1:v$d1,",",col,"]"))
              }
            }
            dd = v$d1 * v$d2
            colnames(RET)[(k-1)*v$dims + 1:dd] <- cnames
          } # end of else of sym
        } # end of for k
      } # end of D==2
      
    } # end of else of isy
  }else{
    # not class-specific
    if(v$isy){
      # not class-specific, but still y-specific
      
      if(v$D==0){
        colnames(RET) <- paste0(p,"_",Ys)
      } # end of D==0
      
      if(v$D==1){
        dimy = 0
        for(y in 1:length(Ys)){
          d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
          colnames(RET)[dimy + 1:d1] <- paste0(p,"_",Ys[y],"[",1:d1,"]")
          dimy = dimy + d1
        }
      } # end of D==1
      
      if(v$D==2){
        dimy = 0
        for(y in 1:length(Ys)){
          d1 = ifelse(missing(yspecd1), v$d1, yspecd1[Ys[y]])
          d2 = ifelse(missing(yspecd2), v$d2, yspecd2[Ys[y]])
          cnames <- c()
          if(v$sym){
            # matrix is symmetrical and only upper-right triangle is stored
            # d1 = d2
            if(v$diag){
              for(col in 1:d1){
                for(row in 1:col){
                  cnames <- c(cnames, paste0(p,"_",Ys[y],"[",row,",",col,"]"))
                }
              }
              dd = d1*(d1+1)/2 # number of elements
              colnames(RET)[dimy + 1:dd] <- cnames
            }else{
              # only non-diagonal elements
              for(col in 2:d1){
                for(row in 1:(col-1)){
                  cnames <- c(cnames, paste0(p,"_",Ys[y],"[",row,",",col,"]"))
                }
              }
              dd = (d1-1)*d1/2 # number of elements
              colnames(RET)[dimy + 1:dd] <- cnames
            } # end of else of diag
          }else{
            # matrix is general rectangular --> all is stored
            if(v$BYROW){
              for(row in 1:d1){
                cnames <- c(cnames, paste0(p,"_",Ys[y],"[",row,",",1:d2,"]"))
              }
            }else{
              for(col in 1:d2){
                cnames <- c(cnames, paste0(p,"_",Ys[y],"[",1:d1,",",col,"]"))
              }
            }
            dd = d1*d2
            colnames(RET)[dimy + 1:(d1*d2)] <- cnames
          } # end of else of sym
          dimy = dimy + dd
        } # end of for y
      } # end of D==2
      
    }else{
      # is NEITHER class-specific, nor y-specific
      
      if(v$D==0){
        colnames(RET) <- p
      } # end of D==0
      
      if(v$D==1){
        colnames(RET) <- paste0(p,"[",1:v$d1,"]")
      } # end of D==1
      
      if(v$D==2){
        cnames <- c()
        if(v$sym){
          # matrix is symmetrical and only upper-right triangle is stored
          # d1 = d2
          if(v$diag){
            for(col in 1:v$d1){
              for(row in 1:col){
                cnames <- c(cnames, paste0(p,"[",row,",",col,"]"))
              }
            }
          }else{
            # only non-diagonal elements
            for(col in 2:v$d1){
              for(row in 1:(col-1)){
                cnames <- c(cnames, paste0(p,"[",row,",",col,"]"))
              }
            }
          } # end of else of diag
        }else{
          # matrix is general rectangular --> all is stored
          if(v$BYROW){
            for(row in 1:v$d1){
              cnames <- c(cnames, paste0(p,"[",row,",",1:v$d2,"]"))
            }
          }else{
            for(col in 1:v$d2){
              cnames <- c(cnames, paste0(p,"[",1:v$d1,",",col,"]"))
            }
          }
        } # end of else of sym
        colnames(RET) <- cnames
      } # end of D==2
    } # end of else of isy
  } # end of else of isspec
  
  return(RET)
}
