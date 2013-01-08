predictEventProb.rsf <- function(object,newdata,times,cause,...){
  ## require(randomSurvivalForest)
  ## class(object) <- c("rsf", "grow")
  ensbCHF <- randomSurvivalForest::predict.rsf(object, test=newdata)
  getCIF <- randomSurvivalForest::competing.risk(ensbCHF, plot=FALSE)$cif.ensb
  if (missing(cause)) cause <- 1
  cif <- getCIF[,,cause]
  Time <- ensbCHF$timeInterest
  pos <- sindex(jump.times=Time,eval.times=times)
  p <- cbind(0,cif)[,pos+1]
  if (is.null(dim(p)))
    {if (length(p)!=length(times))
       stop("Prediction failed")}
  else{
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
  }
  p
}
