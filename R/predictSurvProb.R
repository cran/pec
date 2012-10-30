predictSurvProb <- function(object,newdata,times,...){
  UseMethod("predictSurvProb",object)
}

predictSurvProb.default <- function(object,newdata,times,...){
  stop("No method for evaluating predicted probabilities from objects in class: ",class(object),call.=FALSE)
}


predictSurvProb.numeric <- function(object,newdata,times,...){
  if (NROW(object) != NROW(newdata) || NCOL(object) != length(times))
    stop("Prediction failed")
  object
}

predictSurvProb.matrix <- function(object,newdata,times,...){
  if (NROW(object) != NROW(newdata) || NCOL(object) != length(times)){
    stop(paste("Prediction matrix has wrong dimensions: ",NROW(object)," rows and ",NCOL(object)," columns.\n But requested are predicted probabilities for ",NROW(newdata), " subjects (rows) in newdata and ",NCOL(newdata)," time points (columns)",sep=""))
  }
  object
}

predictSurvProb.aalen <- function(object,newdata,times,...){
  ## require(timereg)
  time.coef <- data.frame(object$cum)
  ntime <- nrow(time.coef)
  objecttime <- time.coef[,1,drop=TRUE]
  ntimevars <- ncol(time.coef)-2
  covanames <- names(time.coef)[-(1:2)]
  notfound <- match(covanames,names(newdata),nomatch=0)==0
  if (any(notfound))
    stop("\nThe following predictor variables:\n\n",
         paste(covanames[notfound],collapse=","),
         "\n\nwere not found in newdata, which only provides the following variables:\n\n",
         paste(names(newdata),collapse=","),
         "\n\n")
  time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
  nobs <- nrow(newdata)
  hazard <- .C("survest_cox_aalen",
               timehazard=double(ntime*nobs),
               as.double(unlist(time.coef[,-1])),
               as.double(unlist(time.vars)),
               as.integer(ntimevars+1),
               as.integer(nobs),
               as.integer(ntime),PACKAGE="pec")$timehazard
  hazard <- matrix(hazard,ncol=ntime,nrow=nobs,dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
  surv <- pmin(exp(-hazard),1)
  if (missing(times)) times <- sort(unique(objecttime))
  pred <- surv[,sindex(jump.times=objecttime,eval.times=times)]
  pred
  if (NROW(pred) != NROW(newdata) || NCOL(pred) != length(times))
    stop("Prediction failed")
  pred
}

predictSurvProb.cox.aalen <- function(object,newdata,times,...){
  #  require(timereg)
  "survest.cox.aalen" <- function(fit,newdata,times,...){
    ##  The time-constant effects first
    const <- c(fit$gamma)
    names(const) <- substr(dimnames(fit$gamma)[[1]],6,nchar(dimnames(fit$gamma)[[1]])-1)
    constant.part <- t(newdata[,names(const)])*const
    constant.part <- exp(colSums(constant.part))
    ##  Then extract the time-varying effects
    time.coef <- data.frame(fit$cum)
    ntime <- nrow(time.coef)
    fittime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    time.part <- .C("survest_cox_aalen",timehazard=double(ntime*nobs),as.double(unlist(time.coef[,-1])),as.double(unlist(time.vars)),as.integer(ntimevars+1),as.integer(nobs),as.integer(ntime),PACKAGE="pec")$timehazard
    time.part <- matrix(time.part,
                        ncol=ntime,
                        nrow=nobs,
                        dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-time.part*constant.part),1)
    if (missing(times)) times <- sort(unique(fittime))
    pred <- surv[,sindex(fittime,times)]
    class(pred) <- c("survest","cox.aalen")
    pred
  }
  p <- survest.cox.aalen(object,times=times,newdata=newdata)
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}


predictSurvProb.mfp <- function(object,newdata,times,...){
  # require(mfp)
  p <- predictSurvProb.coxph(object$fit,newdata=newdata,times=times)
  p
}


predictSurvProb.survnnet <- function(object,newdata,times,train.data,...){
#predictSurvProb.survnnet <- function(object,newdata,times,...){
  require(rms)
  learndat <- train.data
  learndat$nnetFactor <- predict(object,train.data,...)
  newdata$nnetFactor <- predict(object,newdata)
  nnet.form <- reformulate("nnetFactor",object$call$formula[[2]])
  fit.nnet <- cph(nnet.form,data=learndat,se.fit=FALSE,surv=TRUE,x=TRUE,y=TRUE)
  p <- predictSurvProb.cph(fit.nnet,newdata=newdata,times=times)
  p
}


predictSurvProb.rpart <- function(object,newdata,times,train.data,...){
#  require(rpart)
  require(rms)
  learndat <- train.data
  nclass <- length(unique(object$where))
  learndat$rpartFactor <- factor(predict(object,newdata=train.data,...))
  newdata$rpartFactor <- factor(predict(object,newdata=newdata))
  rpart.form <- reformulate("rpartFactor",eval(object$call$formula)[[2]])  
  ##   rpart.form <- reformulate("rpartFactor",object$call$formula[[2]])
  #  fit.rpart <- cph(rpart.form,data=learndat,se.fit=FALSE,surv=TRUE,x=TRUE,y=TRUE)
  fit.rpart <- prodlim(rpart.form,data=learndat)
  p <- predictSurvProb(fit.rpart,newdata=newdata,times=times)
  #  print(p[100:113,1:10])
  p
}


predictSurvProb.coxph <- function(object,newdata,times,...){
  ## baselineHazard.coxph(object,times)
  ## require(survival)
  ## new feature of the survival package requires that the
  ## original data are included
  survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
  survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
  survfit.object <- survival.survfit.coxph(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  inflated.pred <- survival.summary.survfit(survfit.object,times=times)$surv
  p <- t(inflated.pred)
  if ((miss.time <- (length(times) - NCOL(p)))>0)
    p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

predictSurvProb.coxph.penal <- function(object,newdata,times,...){
  ## require(survival)
  frailhistory <- object$history$'frailty(cluster)'$history
  frailVar <- frailhistory[NROW(frailhistory),1]
  ## survfit.object <- survival.survfit.coxph(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  linearPred <- predict(object,newdata=newdata,se.fit=FALSE,conf.int=FALSE)
  basehaz <- basehaz(object)
  bhTimes <- basehaz[,2]
  bhValues <- basehaz[,1]
  survPred <- do.call("rbind",lapply(1:NROW(newdata),function(i){
    (1+frailVar*bhValues*exp(linearPred[i]))^{-1/frailVar}
  }))
  where <- sindex(jump.times=bhTimes,eval.times=times)
  p <- cbind(1,survPred)[,where+1]
  if ((miss.time <- (length(times) - NCOL(p)))>0)
    p <- cbind(p,matrix(rep(NA,miss.time*NROW(p)),nrow=NROW(p)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}



predictSurvProb.cph <- function(object,newdata,times,...){
  if (!match("surv",names(object),nomatch=0)) stop("Argument missing: set surv=TRUE in the call to cph!")
  p <- survest(object,times=times,newdata=newdata,se.fit=FALSE,what="survival")$surv
  if (is.null(dim(p))) p <- matrix(p,nrow=NROW(newdata))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

predictSurvProb.prodlim <- function(object,newdata,times,...){
  require(prodlim)
  p <- predict(object=object,
               type="surv",
               newdata=newdata,
               times=times,
               mode="matrix",
               level.chaos=1)
  if (NROW(newdata)==1 && class(p)=="list"){
    p <- unlist(p)
  }
  if (is.null(dim(p)) && NROW(newdata)>1){
    ## if the model has no covariates
    ## then all cases get the same prediction
    ## in this exceptional case we proceed a vector
    ## p[is.na(p)] <- 0
    p <- as.vector(p)
    if (length(p)!=length(times))
      stop("Prediction failed")
  }
  else{
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
      stop("Prediction failed")
  }
  p
}

predict.survfit <- function(object,newdata,times,bytimes=TRUE,fill="last",...){
    if (length(class(object))!=1 || class(object)!="survfit" || object$typ !="right")
      stop("Predictions only available \nfor class 'survfit', possibly stratified Kaplan-Meier fits.\n For class 'cph' Cox models see survest.cph.")
    
  if (missing(newdata))
    npat <- 1
  else
    if (is.data.frame(newdata))
      npat <- nrow(newdata)
    else stop("If argument `newdata' is supplied it must be a dataframe." )
  
  ntimes <- length(times)
    
  sfit <- summary(object,times=times)

  if (is.na(fill))
    Fill <- function(x,len){x[1:len]}
  else if (fill=="last")
    Fill <- function(x,len){
      y <- x[1:len]
      y[is.na(y)] <- x[length(x)]
      y}
  else stop("Argument fill must be the string 'last' or NA.")

  if (is.null(object$strata)){
    p <- Fill(sfit$surv,ntimes)
    pred <- matrix(rep(p,npat),
                   ncol=ifelse(bytimes,ntimes,npat),
                   nrow=ifelse(bytimes,npat,ntimes),
                   byrow=bytimes)
  }
  else{
    covars <- attr(terms(eval.parent(object$call$formula)),"term.labels")
    ## if (!all(match(covars,names(newdata),nomatch=FALSE)))
      stop("Not all strata defining variables occur in newdata.")

    ## FIXME there are different ways to build strata levels
    ## how can we test which one was used???
    stratdat <- newdata[,covars,drop=FALSE]
    names(stratdat) <- covars
    
    NewStratVerb <- strata(stratdat)
    NewStrat <- interaction(stratdat,sep=" ")
    
    levs <- levels(sfit$strata)
    #    print(levs)
    #    print(levels(NewStrat))
    #    print(levels(NewStratVerb))
    if (!all(choose <- match(NewStratVerb,levs,nomatch=F))
        &&
        !all(choose <- match(NewStrat,levs,nomatch=F)))
      stop("Not all strata levels in newdata occur in fit.")
    survlist <- split(sfit$surv,sfit$strata)
    
    p <- lapply(survlist[choose],Fill,ntimes)
    pred <- matrix(unlist(p,use.names=FALSE),
                   ncol=ifelse(bytimes,ntimes,npat),
                   nrow=ifelse(bytimes,npat,ntimes),
                   byrow=bytimes)
  }
  pred
}

predictSurvProb.survfit <- function(object,newdata,times,...){
  p <- predict.survfit(object,newdata=newdata,times=times,bytimes=TRUE,fill="last")
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}


## library randomSurvivalForest
## predictSurvProb.rsf <- function(object,newdata,times,...){
##   p <- predict.rsf(object,newdata=newdata,times=times,bytimes=TRUE,fill="last")
##   if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
##     stop("Prediction failed")
##   p
## }


predictSurvProb.psm <- function(object,newdata,times,...){
  p <- survest(object,times=times,newdata=newdata,what="survival",conf.int=FALSE)
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

predictSurvProb.phnnet <- function(object,newdata,times,train.data,...){
#  require(survnnet)
  learndat <- train.data
  seeds <- sample(1:1000,size=10)
  object$call$data <- learndat
  re.fitter <- lapply(seeds,function(s){
    set.seed(s)
    refit <- eval(object$call)
    list(learn=predict(refit,learndat),
         val=predict(refit,newdata))
  })
  learndat$nnetFactor <- rowMeans(do.call("cbind",lapply(re.fitter,function(x)x[["learn"]])))
  newdata$nnetFactor <- rowMeans(do.call("cbind",lapply(re.fitter,function(x)x[["val"]])))
  #  a <- predict(object,train.data,...)
  #  b <- predict(object,newdata)
  #  print(cbind(do.call("cbind",lapply(re.fitter,function(x)x[["learn"]]))[1:10,],learndat$nnetFactor[1:10]))
  #  stop()
  #  learndat$nnetFactor <- predict(object,train.data,...)
  #  newdata$nnetFactor <- predict(object,newdata)
  nnet.form <- reformulate("nnetFactor",object$call$formula[[2]])
  ## require(rms)
  fit.nnet <- cph(nnet.form,data=learndat,se.fit=FALSE,surv=TRUE,x=TRUE,y=TRUE)
  p <- predictSurvProb.cph(fit.nnet,newdata=newdata,times=times)
  p
}


predictSurvProb.riskRegression <- function(object,newdata,times,...){
  if (missing(times))stop("Argument times is missing")
  temp <- predict(object,newdata=newdata)
  pos <- sindex(jump.times=temp$time,eval.times=times)
  p <- cbind(1,1-temp$cuminc)[,pos+1,drop=FALSE]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

predictSurvProb.rfsrc <- function(object, newdata, times, ...){
  ptemp <- predict(object,newdata=newdata,...)$survival
  pos <- sindex(jump.times=object$time.interest,eval.times=times)
  p <- cbind(1,ptemp)[,pos+1]
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}


# methods for uncensored regression
# --------------------------------------------------------------------

predictProb <- function(object,newdata,times,...){
  UseMethod("predictProb",object)
}

predictProb.glm <- function(object,newdata,times,...){
  ## no censoring -- only normal family with mu=0 and sd=sd(y)
  N <- NROW(newdata)
  NT <- length(times)
  if (!(unclass(family(object))$family=="gaussian"))
    stop("Currently only gaussian family implemented for glm.")
  betax <- predict(object,newdata=newdata,se.fit=FALSE)
  ##   print(betax[1:10])
  pred.matrix <- matrix(rep(times,N),byrow=TRUE,ncol=NT,nrow=N)
  p <- 1-pnorm(pred.matrix - betax,mean=0,sd=sqrt(var(object$y)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}


predictProb.ols <- function(object,newdata,times,...){
  ## no censoring -- only normal family with mu=0 and sd=sd(y)
  N <- NROW(newdata)
  NT <- length(times)
  if (!(unclass(family(object))$family=="gaussian"))
    stop("Currently only gaussian family implemented.")
  betax <- predict(object,newdata=newdata,type="lp",se.fit=FALSE)
  ##   print(betax[1:10])
  pred.matrix <- matrix(rep(times,N),byrow=TRUE,ncol=NT,nrow=N)
  p <- 1-pnorm(pred.matrix - betax,mean=0,sd=sqrt(var(object$y)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}

predictProb.randomForest <- function(object,newdata,times,...){
  ## no censoring -- only normal family with mu=0 and sd=sd(y)
  N <- NROW(newdata)
  NT <- length(times)
  predMean <- predict(object,newdata=newdata,se.fit=FALSE)
  pred.matrix <- matrix(rep(times,N),byrow=TRUE,ncol=NT,nrow=N)
  p <- 1-pnorm(pred.matrix - predMean,mean=0,sd=sqrt(var(object$y)))
  if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
    stop("Prediction failed")
  p
}





## update.cox <- function(object,tstar,data){
## object$call$data <- data[data$time>tstar,]
## update <- eval(object$call)
## class(update) <- "dynamicCox"
## update
## }
## predictProb.dynamicCox <- function(object,newdata,cutpoints,learn.data,...){
## p <- matrix(1,nrow=NROW(newdata),ncol=length(cutpoints))
## p
## }
