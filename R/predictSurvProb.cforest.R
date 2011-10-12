# CFOREST
# --------------------------------------------------------------------
pecCforest <- function(formula,data,...){
  ## require(party)
  out <- list(forest=cforest(formula,data,...))
  class(out) <- "pecCforest"
  out$call <- match.call()
  out  
}


predictSurvProb.pecCforest <- function (object, newdata, times, ...) {
  survObj <- treeresponse(object$forest,newdata=newdata)
  p <- do.call("rbind",lapply(survObj,function(x){
    predictSurvProb(x,newdata=newdata[1,,drop=FALSE],times=times)
  }))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop("Prediction failed")
  p
}
