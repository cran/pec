# CTREE
# --------------------------------------------------------------------
pecCtree <- function(...){
 out <- list(ctree=ctree(...))
 class(out) <- "pecCtree"
 out$call <- match.call()
 out  
}

predictSurvProb.pecCtree <- function (object, newdata, times, ...) {
  ## require(party)
  N <- NROW(newdata)
  NT <- length(times)
  survObj <- party::treeresponse(object$ctree, newdata=newdata)
  p <- do.call("rbind", lapply(survObj,function(x){
    predictSurvProb(x, newdata=newdata[1,,drop=FALSE], times=times)
  }))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    stop("Prediction failed")
  p
}
