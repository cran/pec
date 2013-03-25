# methods for competing risk regression
# --------------------------------------------------------------------

predictLifeYearsLost <- function(object,newdata,times,cause,...){
  UseMethod("predictLifeYearsLost",object)
}

predictLifeYearsLost.matrix <- function(object,newdata,times,...){
  if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
    stop(paste("Life-years-lost matrix has wrong dimensions: ",
               NROW(object),
               " rows and ",
               NCOL(object),
               " columns.\n But requested are predicted probabilities for ",
               NROW(newdata),
               " subjects (rows) in newdata and ",
               length(times),
               " time points (columns)",
               sep=""))
  }
  object
}

predictLifeYearsLost.prodlim <- function(object,newdata,times,cause,...){
  ## require(prodlim)
  time.interest <- object$time
  cif <- predict(object=object,cause=cause,type="cuminc",newdata=newdata,times=time.interest,mode="matrix",level.chaos=1)
  ## if the model has no covariates
  ## then all cases get the same cif
  ## in this exceptional case we proceed a vector
  if (NROW(cif)==1 && NROW(newdata)>1)
    cif <- as.vector(cif)
  pos <- sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
    stop("Prediction of life-years-lost failed")
  lyl
}

predictLifeYearsLost.FGR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  time.interest <- sort(unique(object$crrFit$uftime))
  cif <- predict(object,newdata=newdata,times=time.interest)
  pos <- sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
    stop("Prediction of life-years-lost failed")
  lyl
}

predictLifeYearsLost.riskRegression <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  time.interest <- object$time
  cif <- predict(object,newdata=newdata,times=time.interest)
  pos <- sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
    stop("Prediction of life-years-lost failed")
  lyl
}

predictLifeYearsLost.ARR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  time.interest <- object$time
  cif <- predict(object,newdata=newdata,times=time.interest)
  pos <- sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
    stop("Prediction of life-years-lost failed")
  lyl
}


predictLifeYearsLost.CauseSpecificCox <- function (object, newdata, times, cause, ...) {
  survtype <- object$survtype
  N <- NROW(newdata)
  NC <- length(object$model)
  eTimes <- object$eventTimes
  if (missing(cause))
    cause <- object$theCause
  causes <- object$causes
  stopifnot(match(as.character(cause),causes,nomatch=0)!=0)
  if (survtype=="survival"){
    if (object$theCause!=cause)
      stop("Object can be used to predict cause ",object$theCause," but not ",cause,".\nNote: the cause can be specified in CSC(...,cause=).")
  }
  # predict cumulative cause specific hazards
  cumHaz1 <- -log(predictSurvProb(object$models[[paste("Cause",cause)]],times=eTimes,newdata=newdata))
  if (length(eTimes)==1)
    Haz1 <- cumHaz1
  else
    Haz1 <- t(apply(cbind(0,cumHaz1),1,diff))
  if (survtype=="hazard"){
    cumHazOther <- lapply(causes[-match(cause,causes)],function(c){
      -log(predictSurvProb(object$models[[paste("Cause",c)]],times=eTimes,newdata=newdata))
    })
    lagsurv <- exp(-cumHaz1 - Reduce("+",cumHazOther))
    cif <- t(apply(lagsurv*Haz1,1,cumsum))
  }
  else{
    tdiff <- min(diff(eTimes))/2
    lagsurv <- pec:::predictSurvProb(object$models[["OverallSurvival"]],times=eTimes-tdiff,newdata=newdata)
    cif <- t(apply(lagsurv*Haz1,1,cumsum))
  }
  pos <- sindex(jump.times=eTimes,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, eTimes)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
    stop("Prediction of life-years-lost failed")
  lyl
}




predictLifeYearsLost.coxboost <- function(object,newdata,times,cause,...){
  if (missing(cause)) stop("missing cause")
  ## if (cause!=1) stop("CoxBoost can only predict cause 1")
  if (attr(object$response,"model")!="competing.risks") stop("Not a competing risk object")
  newcova <- model.matrix(terms(object$formula,data=newdata),
                          data=model.frame(object$formula,data=newdata))[,-c(1)]
  time.interest <- sort(unique(object$coxboost$time))
  cif <- predict(object$coxboost,newdata=newcova,type="CIF",times=time.interest)
  pos <- sindex(jump.times=time.interest,eval.times=times)
  lyl <- matrix(unlist(lapply(1:length(pos), function(j) {
    pos.j <- 1:(pos[j]+1)
    p <- cbind(0,cif)[,pos.j,drop=FALSE]
    time.diff <- diff(c(0, object$time.interest)[pos.j])
    apply(p, 1, function(x) {sum(x[-length(x)] * time.diff)})
  })), ncol = length(pos))
  if (NROW(lyl) != NROW(newdata) || NCOL(lyl) != length(times))
    stop("Prediction of life-years-lost failed")
  lyl
}
