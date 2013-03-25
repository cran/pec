ipcw <- function(formula,
                 data,
                 method=c("cox","marginal","nonpar","aalen","none"),
                 times,
                 subjectTimes,
                 subjectTimesLag=1,
                 what){
  if (!missing(what))
    stopifnot(all(match(what,c("IPCW.times","IPCW.subjectTimes"))))
  method <- match.arg(method,c("cox","marginal","nonpar","aalen","none"))
  class(method) <- method
  UseMethod("ipcw",method)
}


ipcw.none <- function(formula,data,method,times,subjectTimes,subjectTimesLag,what){
  if (missing(subjectTimesLag)) subjectTimesLag=1
  if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
  call <- match.call()
  environment(call$formula) <- NULL
  #  weigths at requested times
  if (match("IPCW.times",what,nomatch=FALSE)){
    length.times <- length(times)
    stopifnot(length.times>0)
    IPCW.times <- rep(1,length(times))
    if (is.null(dim(IPCW.times))){
      names(IPCW.times) <- NULL
    }
    else{
      names(IPCW.times) <- paste("t",1:length(IPCW.times),sep="")
    }
  }
  else
    IPCW.times <- NULL
  #  weigths at subject specific event times
  if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
    IPCW.subjectTimes <- rep(1,length(subjectTimes))
    names(IPCW.subjectTimes) <- paste("T",1:length(IPCW.subjectTimes),".lag",subjectTimesLag,sep="")
  }
  else
    IPCW.subjectTimes <- NULL
  fit <- NULL
  out <- list(times=times,
              IPCW.times=IPCW.times,
              IPCW.subjectTimes=IPCW.subjectTimes,
              fit=fit,
              call=call,
              method=method)
  class(out) <- "IPCW"
  out
}
## reverse Random Survival Forests
## ipcw.rfsrc <- function(formula,data,method,times,subjectTimes,subjectTimesLag,what){
  ## if (missing(subjectTimesLag)) subjectTimesLag=1
  ## if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
  ## call <- match.call()
  ## require(randomForestSRC)
  ## status.name <- all.vars(formula)[2]
  ## reverse.data <- data
  ## reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  ## stopifnot(NROW(na.omit(data))>0)
  ## fit <- rsfrc(formula,data=reverse.data,surv=TRUE,x=TRUE,y=TRUE)
  ## #  weigths at requested times
  ## if (match("IPCW.times",what,nomatch=FALSE)){
    ## length.times <- length(times)
    ## stopifnot(length.times>0)
    ## IPCW.times <- predictSurvProb(fit,newdata=data,times=times)
    ## if (is.null(dim(IPCW.times))){
      ## names(IPCW.times) <- NULL
    ## } else{
      ## colnames(IPCW.times) <- paste("t",1:NCOL(IPCW.times),sep="")
      ## rownames(IPCW.times) <- 1:NROW(IPCW.times)
    ## }
  ## }
  ## else
    ## IPCW.times <- NULL
  ## #  weigths at subject specific event times
  ## if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
    ## if (subjectTimesLag==1)
      ## IPCW.subjectTimes <- predictSurvProb(fit,times=subjectTimes-min(diff(c(0,unique(subjectTimes))))/2,newdata=data)
    ## else if (subjectTimesLag==0){
      ## IPCW.subjectTimes <- survest(fit,times=subjectTimes,what='parallel')
    ## }
    ## else stop("SubjectTimesLag must be 0 or 1")
    ## names(IPCW.subjectTimes) <- paste("T",1:length(IPCW.subjectTimes),".lag",subjectTimesLag,sep="")
  ## }
  ## else
    ## IPCW.subjectTimes <- NULL
  ## out <- list(times=times,
              ## IPCW.times=IPCW.times,
              ## IPCW.subjectTimes=IPCW.subjectTimes,
              ## fit=fit,
              ## call=call,
              ## method=method)
  ## class(out) <- "IPCW"
  ## out
## }


## reverse Kaplan-Meier 
ipcw.marginal <- function(formula,data,method,times,subjectTimes,subjectTimesLag,what){
  if (missing(subjectTimesLag)) subjectTimesLag=1
  if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
  call <- match.call()
  environment(call$formula) <- NULL
  formula <- update.formula(formula,"~1")
  fit <- prodlim(formula,data=data,reverse=TRUE)
  #  weigths at requested times
  if (match("IPCW.times",what,nomatch=FALSE)){
    length.times <- length(times)
    stopifnot(length.times>0)
    IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
    if (length(times)==1){
      names(IPCW.times) <- NULL
    }
    else{
      names(IPCW.times) <- paste("t",1:length(IPCW.times),sep="")
    }
  }
  else
    IPCW.times <- NULL
  #  weigths at subject specific event times
  if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
    IPCW.subjectTimes <- prodlim:::predictSurvIndividual(fit,lag=subjectTimesLag)
    names(IPCW.subjectTimes) <- paste("T",1:length(IPCW.subjectTimes),".lag",subjectTimesLag,sep="")
  }
  else
    IPCW.subjectTimes <- NULL
  out <- list(times=times,IPCW.times=IPCW.times,IPCW.subjectTimes=IPCW.subjectTimes,fit=fit,call=call,method=method)
  class(out) <- "IPCW"
  out
  ##   locsubjectTimes <- match(subjectTimes,fit$time,nomatch=NA)
  ##   if (any(is.na(locsubjectTimes))) stop("Can not locate all individual observation times" )
  ##   IPCW.subjectTimes <- c(1,fit$surv)[locsubjectTimes] ## at (subjectTimes_i-)
  ##   IPCW.times <- c(1,fit$surv)[sindex(jump.times=fit$time,eval.times=times) +1] ## at all requested times
}

## reverse Stone-Beran 
ipcw.nonpar <- function(formula,data,method,times,subjectTimes,subjectTimesLag,what){
  if (missing(subjectTimesLag)) subjectTimesLag=1
  if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
  call <- match.call()
  environment(call$formula) <- NULL
  fit <- prodlim(formula,data=data,reverse=TRUE,bandwidth="smooth")
  #  weigths at requested times
  if (match("IPCW.times",what,nomatch=FALSE)){
    length.times <- length(times)
    stopifnot(length.times>0)
    IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
    if (is.null(dim(IPCW.times))){
      names(IPCW.times) <- NULL
    }
    else{
      colnames(IPCW.times) <- paste("t",1:NCOL(IPCW.times),sep="")
      rownames(IPCW.times) <- 1:NROW(IPCW.times)
    }
  }
  else
    IPCW.times <- NULL
  #  weigths at subject specific event times
  if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
    IPCW.subjectTimes <- predictSurvIndividual(fit,lag=subjectTimesLag)
    names(IPCW.subjectTimes) <- paste("T",1:length(IPCW.subjectTimes),".lag",subjectTimesLag,sep="")
  }
  else
    IPCW.subjectTimes <- NULL
  out <- list(times=times,
              IPCW.times=IPCW.times,
              IPCW.subjectTimes=IPCW.subjectTimes,
              fit=fit,
              call=call,
              method=method)
  class(out) <- "IPCW"
  out
}

#reverse Cox via the survival package
#ipcw.coxph <- function(formula,data,method,times,subjectTimes,subjectTimesLag,what){
#  if (missing(subjectTimesLag)) subjectTimesLag=1
#  if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
#  call <- match.call()
#  length.times <- length(times)
#  stopifnot(length.times>0)
#  require(survival)
#  survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
#  survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
#  status.name <- all.vars(formula)[2]
#  reverse.data <- data
#  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
#  stopifnot(NROW(na.omit(data))>0)
#  fit <- coxph(formula,data=reverse.data)
#  survfit.object <- survival.survfit.coxph(fit,newdata=data,se.fit=FALSE,conf.int=FALSE)
#  if (match("IPCW.times",what,nomatch=FALSE))
#    IPCW.times <- survest(fit,newdata=data,times=times,se.fit=FALSE)$surv
#  else
#    IPCW.times <- NULL
#  if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
#    if (subjectTimesLag==1) 
#      IPCW.subjectTimes <- survest(fit,times=subjectTimes-min(diff(c(0,unique(subjectTimes))))/2,what='parallel')
#    else if (subjectTimesLag==0)
#      IPCW.subjectTimes <- survest(fit,times=subjectTimes,what='parallel')
#    else stop("SubjectTimesLag must be 0 or 1")}
#  else
#    IPCW.subjectTimes <- NULL
#  out <- list(times=times,IPCW.times=IPCW.times,IPCW.subjectTimes=IPCW.subjectTimes,fit=fit,call=call,method=method)
#  class(out) <- "IPCW"
#  out
#}

#reverse Cox via Harrel's package
ipcw.cox <- function(formula,data,method,times,subjectTimes,subjectTimesLag,what){
  if (missing(subjectTimesLag)) subjectTimesLag=1
  if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
  call <- match.call()
  environment(call$formula) <- NULL
  ## require(rms)
  status.name <- all.vars(formula)[2]
  reverse.data <- data
  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  stopifnot(NROW(na.omit(data))>0)
  fit <- cph(formula,data=reverse.data,surv=TRUE,x=TRUE,y=TRUE)
  #  weigths at requested times
  if (match("IPCW.times",what,nomatch=FALSE)){
    length.times <- length(times)
    stopifnot(length.times>0)
    IPCW.times <- survest(fit,newdata=data,times=times,se.fit=FALSE)$surv
    if (is.null(dim(IPCW.times))){
      names(IPCW.times) <- NULL
    } else{
      colnames(IPCW.times) <- paste("t",1:NCOL(IPCW.times),sep="")
      rownames(IPCW.times) <- 1:NROW(IPCW.times)
    }
  }
  else
    IPCW.times <- NULL
  #  weigths at subject specific event times
  if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
    if (subjectTimesLag==1)
      IPCW.subjectTimes <- survest(fit,times=subjectTimes-min(diff(c(0,unique(subjectTimes))))/2,what='parallel')
    else if (subjectTimesLag==0){
      IPCW.subjectTimes <- survest(fit,times=subjectTimes,what='parallel')
    }
    else stop("SubjectTimesLag must be 0 or 1")
    names(IPCW.subjectTimes) <- paste("T",1:length(IPCW.subjectTimes),".lag",subjectTimesLag,sep="")
  }
  else
    IPCW.subjectTimes <- NULL
  out <- list(times=times,
              IPCW.times=IPCW.times,
              IPCW.subjectTimes=IPCW.subjectTimes,
              fit=fit,
              call=call,
              method=method)
  class(out) <- "IPCW"
  out
}

#reverse Aalen method via the timereg package
ipcw.aalen <- function(formula,data,method,times,subjectTimes,subjectTimesLag,what){
  if (missing(subjectTimesLag)) subjectTimesLag=1
  if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
  call <- match.call()
  environment(call$formula) <- NULL
  require(timereg)
  ## require(rms)
  status.name <- all.vars(formula)[2]
  reverse.data <- data
  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  fit <- do.call(method,list(formula=formula,data=reverse.data,n.sim=0))
  #  weigths at requested times
  if (match("IPCW.times",what,nomatch=FALSE)){
    length.times <- length(times)
    stopifnot(length.times>0)
    IPCW.times <- predictSurvProb(fit,newdata=data,times=times)
    if (is.null(dim(IPCW.times))){
      names(IPCW.times) <- NULL
    }
    else{
      colnames(IPCW.times) <- paste("t",1:NCOL(IPCW.times),sep="")
      rownames(IPCW.times) <- 1:NROW(IPCW.times)
    }
  }
  else
    IPCW.times <- NULL
  if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
    if (subjectTimesLag==1) 
      IPCW.subjectTimes <- diag(predictSurvProb(fit,newdata=data,times=pmax(0,subjectTimes-min(diff(unique(subjectTimes)))/2)))
    else if (subjectTimesLag==0)
      IPCW.subjectTimes <- diag(predictSurvProb(fit,newdata=data,times=subjectTimes))
    else stop("SubjectTimesLag must be 0 or 1")
    names(IPCW.subjectTimes) <- paste("T",1:length(IPCW.subjectTimes),".lag",subjectTimesLag,sep="")
  }
  else
    IPCW.subjectTimes <- NULL
  out <- list(times=times,IPCW.times=IPCW.times,IPCW.subjectTimes=IPCW.subjectTimes,fit=fit,call=call,method=method)
  class(out) <- "IPCW"
  out
}

#Stone-Beran using the linear predictor of a reverse Cox method
#ipcw.project <- function(formula,data,method,times,subjectTimes){ 
#  if (missing(subjectTimesLag)) subjectTimesLag=1
#  if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
#  call <- match.call()
#  length.times <- length(times)
#  stopifnot(length.times>0)
#  require(rms)
#  status.name <- all.vars(formula)[2]
#  reverse.data <- data
#  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
#  stopifnot(NROW(na.omit(data))>0)
#  fit <- cph(formula,data=reverse.data,surv=TRUE,x=TRUE,y=TRUE)
#  data$lp <- predict(fit,type="lp")
#  reform <- reformulate("lp",formula[[2]])
#  fit <- prodlim(reform,data=data,reverse=TRUE,bandwidth="smooth")
#  IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
#  IPCW.subjectTimes <- predictSurvIndividual(fit,subjectTimesLag=1)
#  out <- list(times=times,IPCW.times=IPCW.times,IPCW.subjectTimes=IPCW.subjectTimes,fit=fit,call=call,method=method)
#  class(out) <- "IPCW"
#  out
#}
