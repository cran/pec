internalReevalFit <- function(object,data,step,silent=FALSE,verbose=FALSE){
  object$call$data <- data
  try2fit <- try(refit <- eval(object$call),silent=silent)
  if (inherits(try2fit,"try-error")){
    if (verbose==TRUE)
      warning(paste("During bootstrapping: model ",class(object)[[1]]," failed in step ",step),immediate.=TRUE)
    NULL
  }
  else
    refit
}
