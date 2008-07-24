"print.pec" <- function(x,
                        times,
                        digits=3,
                        what=NULL,
                        ...){
  
  cat("\nPrediction error curves\n\n")
  prmatrix(matrix(names(x$models),ncol=1,dimnames=list(c(1:length(names(x$models))),"Prediction models")),quote=FALSE)
  if (!is.null(x$cens.model)){
    if (x$cens.model!="none")
      cat("\nIPCW:",x$cens.model,"model")
    else cat("\nno censoring")}
  if (!is.null(x$method)) print(x$method)
  if (missing(times)) times <- x$maxtime ## times <- quantile(x$time,.9)
  
  cat("\nCumulative prediction error between",x$start,"and",times,":\n\n")
  out <- crps(object=x,times=times,start=x$start,what=what)
  if (is.matrix(out))
    print(out,digits=digits,quote=FALSE)
  else{
    print.listof(out,digits=digits,quote=FALSE)
  }
  failed <- sapply(x$failed,length)>0
  if (sum(failed)>0){
    cat("\nWarning: In some bootstrap samples the model failed:\n")
    print(sapply(x$failed[failed],table))
  }
  invisible(x)
}
