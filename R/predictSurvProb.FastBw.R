selectCox <- function(formula,data,rule="aic"){
  ## require(rms)
  ## require(prodlim)
  fit <- cph(formula, data, surv=TRUE)
  bwfit <- fastbw(fit,rule=rule)
  if (length(bwfit$names.kept)==0){
    newform <- reformulate("1",formula[[2]])
    newfit <- prodlim(newform,data=data)
  }
  else{
    newform <- reformulate(bwfit$names.kept, formula[[2]])
    newfit <- cph(newform,data, surv=TRUE)
  }
  out <- list(fit=newfit,In=bwfit$names.kept)
  out$call <- match.call()
  class(out) <- "selectCox"
  out
}

predictSurvProb.selectCox <- function(object,newdata,times,...){
   predictSurvProb(object[[1]],newdata=newdata,times=times,...)
 }
