
# methods for competing risk regression
# --------------------------------------------------------------------

predictEventProb <- function(object,newdata,times,cause,...){
  UseMethod("predictEventProb",object)
}

predictEventProb.matrix <- function(object,newdata,times,...){
  if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
    stop(paste("Prediction matrix has wrong dimensions: ",NROW(object)," rows and ",NCOL(object)," columns.\n But requested are predicted probabilities for ",NROW(newdata), " subjects (rows) in newdata and ",NCOL(newdata)," time points (columns)",sep=""))
  }
  object
}

predictEventProb.prodlim <- function(object,newdata,times,cause,...){
  require(prodlim)
  p <- predict(object=object,cause=cause,type="cuminc",newdata=newdata,times=times,mode="matrix",level.chaos=1)
  if (NROW(p)==1) p <- as.vector(p)
  p[is.na(p)] <- 0
  if (is.null(dim(p)))
    {if (length(p)!=length(times))
       stop("Prediction failed")}
  else{
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
  }
  p
}

predictEventProb.FGR <- function(object,newdata,times,cause,...){
  ## require(cmprsk)
  # require(compRisk)
  p <- predict(object=object,newdata=newdata,times=times)
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

predictEventProb.riskRegression <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  p <- temp$risk
  pos <- sindex(jump.times=temp$time,eval.times=times)
  cbind(0,p)[,pos+1,drop=FALSE]
}

predictEventProb.ARR <- function(object,newdata,times,cause,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata,times=times)
  p <- temp$P1
  pos <- sindex(jump.times=temp$time,eval.times=times)
  cbind(0,p)[,pos+1,drop=FALSE]
}


predictEventProb.CauseSpecificCox <- function (object, newdata, times, cause, ...) {
  survtype <- object$survtype
  N <- NROW(newdata)
  NC <- length(object$model)
  eTimes <- object$eventTimes
  if (missing(cause))
    cause <- object$theCause
  causes <- object$causes
  stopifnot(match(as.character(cause),causes,nomatch=0)!=0)
  # predict cumulative cause specific hazards
  cumHaz1 <- -log(predictSurvProb(object$models[[paste("Cause",cause)]],
                                  times=eTimes,
                                  newdata=newdata))
  if (length(eTimes)==1)
    Haz1 <- cumHaz1
  else
    Haz1 <- t(apply(cbind(0,cumHaz1),1,diff))
  if (survtype=="hazard"){
    cumHazOther <- lapply(causes[-match(cause,causes)],function(c){
      -log(predictSurvProb(object$models[[paste("Cause",c)]],times=eTimes,newdata=newdata))
    })
    lagsurv <- exp(-cumHaz1- do.call("+",cumHazOther))
    cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  }
  else{
    tdiff <- min(diff(eTimes))/2
    lagsurv <- pec:::predictSurvProb(object$models[["OverallSurvival"]],times=eTimes-tdiff,newdata=newdata)
    cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  }
  pos <- sindex(jump.times=eTimes, eval.times=times)
  cbind(0,cuminc1)[,pos+1,drop=FALSE]
}

## predictUpdateProb.CSC <- function (object, newdata,times,horizon, cause, ...) {
  ## survtype <- object$survtype
  ## N <- NROW(newdata)
  ## NC <- length(object$model)
  ## eTimes <- object$eventTimes
  ## if (missing(cause))
    ## cause <- object$theCause
  ## causes <- object$causes
  ## stopifnot(match(as.character(cause),causes,nomatch=0)!=0)
  ## # predict cumulative cause specific hazards
  ## cumHaz1 <- -log(predictSurvProb(object$models[[paste("Cause",cause)]],times=eTimes,newdata=newdata))
  ## Haz1 <- t(apply(cbind(0,cumHaz1),1,diff))
  ## if (survtype=="hazard"){
    ## cumHazOther <- lapply(causes[-match(cause,causes)],function(c){
      ## -log(predictSurvProb(object$models[[paste("Cause",c)]],times=eTimes,newdata=newdata))
    ## })
    ## lagsurv <- exp(-cumHaz1- do.call("+",cumHazOther))
    ## cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  ## }
  ## else{
    ## tdiff <- min(diff(eTimes))/2
    ## lagsurv <- pec:::predictSurvProb(object$models[["OverallSurvival"]],times=eTimes-tdiff,newdata=newdata)
    ## cuminc1 <- t(apply(lagsurv*Haz1,1,cumsum))
  ## }
  ## pos <- sindex(jump.times=eTimes, eval.times=times)
  ## cbind(0,cuminc1)[,pos+1,drop=FALSE]
## }
