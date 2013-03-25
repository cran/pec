## rsf.default=randomSurvivalForest:::rsf.default 
## predictSurvProb.rsf <- function (object, newdata, times, ...)  { 
  ## N <- NROW(newdata) 
  ## ## class(object) <- c("rsf", "grow")
  ## S <- exp(-predict.rsf(object, test=newdata)$ensemble)
  ## if (N==1) S <- matrix(S,nrow=1)
  ## Time <- object$timeInterest 
  ## p <- cbind(1,S)[,1+sindex(Time, times),drop=FALSE] 
  ## if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))  
    ## stop("Prediction failed") 
  ## p 
## } 



