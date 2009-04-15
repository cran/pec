"print.Cindex" <- function(x,digits=3,what,...){
  cat("\nThe c-index quantifies the ability
 of a model to order the risk of individuals
 according to the individual event times. \n\n")
  prmatrix(matrix(names(x$models),ncol=1,dimnames=list(c(1:length(names(x$models))),"Prediction models")),quote=FALSE)
  if (!is.null(x$method)) print(x$method)
  res <- c("PredCindex","AppCindex","OutOfBagCindex","NoInfCindex")
  names(res) <- c(x$method$name,c("App","Oob","NoInf"))
  found <- match(names(x),res,nomatch=FALSE)
  res <- res[found]
  if (!missing(what)){
    res <- res[names(res) %in% what]
  }
  out <- lapply(res,function(r){
    out <- sapply(1:length(x$models),function(w){
      x[[r]][[w]]
    })
    names(out) <- names(x$models)
    out
  })
  outMat <- 100*do.call("cbind",out)
  outMat <- outMat[order(outMat[,NCOL(outMat)]),,drop=FALSE]
  cat("\n\nEstimated C-index in %\n\n")
  print(outMat,digits)
  if(x$method$name=="BootCV")
    cat("\n\nApp    : Apparent performance\nBootCV : Bootstrap crossvalidated performance\n\n")
  invisible(out)
}
