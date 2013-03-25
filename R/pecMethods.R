pec <- function(object,...){
  UseMethod("pec",object=object)
}

# {{{ header pec.list

pec.list <- function(object,
                     formula,
                     data,
                     times,
                     cause,
                     start,
                     maxtime,
                     exact=TRUE,
                     exactness=100,
                     fillChar=NA,
                     cens.model="cox",
                     ipcw.refit=FALSE,
                     splitMethod="none",
                     B,
                     M,
                     reference=TRUE,
                     model.args=NULL,
                     model.parms=NULL,
                     keep.index=FALSE,
                     keep.matrix=FALSE,
                     keep.models="Call",
                     keep.residuals=FALSE,
                     keep.pvalues=FALSE,
                     noinf.permute=FALSE,
                     multiSplitTest=FALSE,
                     testIBS,
                     testTimes,
                     confInt=FALSE,
                     confLevel=0.95,
                     verbose=TRUE,
                     savePath=NULL,
                     slaveseed=NULL,
                     ...)
{

  # }}}
  # {{{ checking integrity some arguments

  theCall=match.call()
  if (match("replan",names(theCall),nomatch=FALSE))
    stop("Argument name 'replan' has been replaced by 'splitMethod'.")
  
  if (!missing(testIBS) && (!(is.logical(testIBS) || (length(testIBS)==2 && is.numeric(testIBS)))))
    stop("Argument testIBS can be TRUE/FALSE or a vector of two numeric values.")
  if (missing(testIBS)) testIBS <- FALSE
  if (keep.residuals && missing(testTimes))
    stop("To keep.residuals please specify testTimes.")
  if (missing(splitMethod) && multiSplitTest==TRUE){
    stop("Need data splitting to compute van de Wiel's test")
  }
  if (missing(M) && multiSplitTest) M <- NA

  # }}}
  # {{{ formula

  if (missing(formula)){
    ftry <- try(formula <- eval(object[[1]]$call$formula),silent=TRUE)
    if ((class(ftry)=="try-error") || match("formula",class(formula),nomatch=0)==0)
      stop("Argument formula is missing.")
    else if (verbose)
      warning("Formula missing. Using formula from first model")
  }
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[1]=="~")
      ||
      (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
    stop("Invalid specification of formula.\n Could be that you forgot the right hand side:\n ~covariate1 + covariate2 + ...?\nNote that any subsetting, ie data$var or data[,\"var\"], is not supported by this function.")
  }
  else{
    if (!(formula.names[2] %in% c("Surv","Hist")))
      survp <- FALSE
    else
      survp <- TRUE
  }

  # }}}
  # {{{ data
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (match("data.frame",class(data),nomatch=0)==0)
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I use the data from the call to the first model instead.")
  }
  # }}}
  # {{{ censoring model
  
  cens.model <- match.arg(cens.model,c("cox","marginal","nonpar","aalen","none"))
  
  # }}}
  # {{{ response

  histformula <- formula
  if (histformula[[2]][[1]]==as.name("Surv")){
    histformula[[2]][[1]] <- as.name("Hist")
  }
  m <- model.frame(histformula,data,na.action=na.fail)
  response <- model.response(m)
  if (match("Surv",class(response),nomatch=0)!=0){
    attr(response,"model") <- "survival"
    attr(response,"cens.type") <- "rightCensored"
    model.type <- "survival"
  }
  model.type <- attr(response,"model")
  if (model.type=="competing.risks"){
    predictHandlerFun <- "predictEventProb"
    if (missing(cause))
      cause <- attr(response,"state")[1]
  }
  else{
    if (survp==FALSE && NCOL(response)!=1) stop("Response must be one-dimensional.")
    if (survp==TRUE && NCOL(response)!=2) stop("Survival response must have two columns: time and status.")
    predictHandlerFun <- "predictSurvProb"
  }

  # }}}
  # {{{ prediction models
  if (reference==TRUE) {
    ProdLimform <- as.formula(update(formula,".~NULL"),env=NULL)
    ProdLimfit <- do.call("prodlim",list(formula=ProdLimform,data=data))
    ProdLimfit$call$data <- NULL
    ProdLimfit$formula <- NULL
    ProdLimfit$call$formula=ProdLimform
    ProdLimfit$formula <- as.formula(ProdLimfit$formula,env=NULL)
    ## print(environment(ProdLimfit$formula))
    if (model.type=="competing.risks")
      object <- c(list(Reference=ProdLimfit),object)
    else
      object <- c(list(Reference=ProdLimfit),object)
  }
  if (is.null(names(object))){
    names(object) <- sapply(object,function(o)class(o)[1])
  }
  else{
    names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
  }
  names(object) <- make.names(names(object),unique=TRUE)
  NF <- length(object) 
  # }}}  
  # {{{ sort the data 

  if (survp){
    neworder <- order(response[,"time"],-response[,"status"])
    if (predictHandlerFun=="predictEventProb"){
      event <- getEvent(response,mode="character")
      event <- event[neworder]
    }
    response <- response[neworder,,drop=FALSE]
    Y <- response[,"time"]
    status <- response[,"status"]
  }
  else{
    cens.model <- "none"
    neworder <- order(response)
    Y <- response[neworder]
    status <- rep(1,length(Y))
  }
  ## for competing risks find the cause of interest.
  if (predictHandlerFun=="predictEventProb"){
    availableCauses <- unique(event)
    if (!match(cause, availableCauses,nomatch=FALSE))
      stop("Cause ",cause," is not among the available causes: ",paste(availableCauses,collapse=", "))
    event <- event==cause
  }
  ##   else{
  ##     event <- NULL
  ##   }
  data <- data[neworder,]
  unique.Y <- unique(Y)
  N <- length(Y)
  NU <- length(unique.Y)
  # }}}
  # {{{ splitMethod

  splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
  B <- splitMethod$B
  ResampleIndex <- splitMethod$index
  k <- splitMethod$k
  do.resample <- !(is.null(ResampleIndex))
  if (keep.matrix==TRUE & !do.resample){
    warning("Argument keep.matrix set to FALSE, since no resampling/crossvalidation is requested.")
    keep.matrix <- FALSE
  }

  # }}}      
  # {{{ find maxtime, start, and jumptimes in the range of the response 
  if (missing(maxtime) || is.null(maxtime))
    maxtime <- unique.Y[NU]
  if (missing(start))
    if (survp==TRUE)
      start <- 0  ## survival times are positive
    else
      start <- min(unique.Y) 
  if (missing(times)){
    if (exact==TRUE)
      times <- unique(c(start,unique.Y))
    else
      times <- seq(start,maxtime,(maxtime - start)/exactness)
  }
  else{
    if (exact==TRUE) 
      times <- sort(c(start,unique(times),unique.Y))
    else
      times <- sort(unique(c(start,times)))
  }
  times <- times[times<=maxtime]
  NT <-  length(times)

  # }}}
  # {{{ IPCW (all equal to 1 without censoring) 
  if((cens.model %in% c("aalen","cox","nonpar"))){
    if (all(as.numeric(status)==1) || sum(status)==N){
      if (verbose)
      message("No censored observations: cens.model coerced to \"none\".")
      cens.model <- "none"
    }
    if (length(attr(terms(formula),"factors"))==0){
      if (verbose==TRUE)
        message("No covariates  specified: Kaplan-Meier for censoring times used for weighting.")
      cens.model <- "marginal"}
  }
  if (predictHandlerFun=="predictEventProb"){
    iFormula <- as.formula(paste("Surv(itime,istatus)","~",as.character(formula)[[3]]))
    iData <- data
    iData$itime <- response[,"time"]
    iData$istatus <- response[,"status"]
    if (ipcw.refit==TRUE)
      stop("pec: internal refitting of censoring distribution not (not yet) supported for competing risks")
    ipcw.call <- NULL
    ipcw <- ipcw(formula=iFormula,data=iData,method=cens.model,times=times,subjectTimes=Y,subjectTimesLag=1)
    ipcw$dim <- if (cens.model %in% c("marginal","none")) 0 else 1
  }
  else{
    if (ipcw.refit==TRUE && splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632"))
      ipcw.call <- list(formula=formula,data=NULL,method=cens.model,times=times,subjectTimes=NULL,subjectTimesLag=1)
    else
      ipcw.call <- NULL
    ipcw <- ipcw(formula=formula,data=data,method=cens.model,times=times,subjectTimes=Y,subjectTimesLag=1)
    ipcw$dim <- if (cens.model %in% c("marginal","none")) 0 else 1
  }
  ## force ipc weights not to exaggerate
  ## weights should not be greater than 1/(sample size)
  ## if (ipcw$dim==1){
  ## ipcw$IPCW.times <- apply(ipcw$IPCW.times,1,function(x)pmax(x,1/N))
  ## }
  ## else{
  ## ipcw$IPCW.times <- pmax(ipcw$IPCW.times,1/N)
  ## }
  ## ipcw$IPCW.subjectTimes <- pmax(ipcw$IPCW.subjectTimes,1/N)
  ## browser()
  
  #  wt <- ipcw$IPCW.times
  #  wt.obs <- ipcw$IPCW.subjectTimes
  #  if (NCOL(wt)>1) {stopifnot(length(wt)==(N*NT))}  else{stopifnot(length(wt)==NT)}

  # }}}
  # {{{ checking the models for compatibility with resampling
  if (do.resample){
    cm <- checkModels(object=object,model.args=model.args,model.parms=model.parms,splitMethod=splitMethod$internal.name)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }
  # }}}
  # {{{ ---------------------------Apparent error---------------------------

  AppErr <- lapply(1:NF,function(f){
    ## message(f)
    fit <- object[[f]]
    extraArgs <- model.args[[f]]
    if (predictHandlerFun=="predictEventProb"){
      pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times,train.data=data,cause=cause),extraArgs))
      if (class(object[[f]])[[1]]=="matrix") pred <- pred[neworder,,drop=FALSE]
      .C("pecCR",pec=double(NT),as.double(Y),as.double(status),as.double(event),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
    }
    else{
      pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times,train.data=data),extraArgs))
      if (class(object[[f]])[[1]]=="matrix") pred <- pred[neworder,,drop=FALSE]
      .C("pec",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
    }
  })

  names(AppErr) <- names(object)

  # }}}
  # {{{------------------------No information error------------------------

  if (splitMethod$internal.name %in% c("Boot632plus")){
    if (verbose==TRUE){
      message("Computing noinformation error using all permutations")
    }
    if (noinf.permute==FALSE){
      NoInfErr <- lapply(1:NF,function(f){
        fit <- object[[f]]
        extraArgs <- model.args[[f]]
        pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=times,train.data=data),extraArgs))
        ## browser()
        extraArgs <- model.args[[f]]
        if (predictHandlerFun=="predictEventProb")
          .C("pec_noinfCR",pec=double(NT),as.double(Y),as.double(status),as.double(event),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
        else
          .C("pec_noinf",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred),as.double(ipcw$IPCW.times),as.double(ipcw$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred))),NAOK=TRUE,PACKAGE="pec")$pec
      })
      names(NoInfErr) <- names(object)
    }else{
      if (verbose==TRUE){
        message("Noinformation error simulation loop (B=",B,")")
      }
      NoInfErrList <- lapply(1:B,function(b){
        if (verbose==TRUE){
          internalTalk(b,B,sign=".")
        }
        responseNames <- colnames(response)
        noinf.b <- data[sample(1:NROW(data),replace=FALSE),-match(responseNames,names(data))]
        noinf.b[,responseNames] <- response
        ipcw.b <- ipcw(formula=formula,data=noinf.b,method=cens.model,times=times,subjectTimes=Y,subjectTimesLag=1)
        noinfPredErr <- lapply(1:NF,function(f){
          fit.b <- internalReevalFit(object=object[[f]],data=noinf.b,step=b,silent=FALSE,verbose=verbose)
          ## fit.b$call <- object[[f]]$call
          extraArgs <- model.args[[f]]

          pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=noinf.b,times=times,train.data=data),extraArgs))
          if (predictHandlerFun=="predictEventProb"){
            pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=noinf.b,times=times,train.data=data,cause=cause),extraArgs))
            .C("pecCR",pec=double(NT),as.double(Y),as.double(status),as.double(event),as.double(times),as.double(pred.b),as.double(ipcw.b$IPCW.times),as.double(ipcw.b$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pec
          }
          else{
            pred.b <- do.call(predictHandlerFun,c(list(object=fit.b,newdata=noinf.b,times=times,train.data=data),extraArgs))
            .C("pec",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred.b),as.double(ipcw.b$IPCW.times),as.double(ipcw.b$IPCW.subjectTimes),as.integer(N),as.integer(NT),as.integer(ipcw$dim),as.integer(is.null(dim(pred.b))),NAOK=TRUE,PACKAGE="pec")$pec
          }
        })
        noinfPredErr
      })
      NoInfErrMat <- lapply(1:NF,function(f){
        do.call("rbind",lapply(NoInfErrList,function(x){
          x[[f]]
        }))})
      NoInfErr <- lapply(NoInfErrMat,colMeans)
      names(NoInfErr) <- names(object)
    }
  }

  # }}}
  # {{{--------------k-fold and leave-one-out CrossValidation-----------------------

  if (splitMethod$internal.name %in% c("crossval","loocv")){
    kCV <- kFoldCrossValidation(object=object,data=data,Y=Y,status=status,event=event,times=times,cause=cause,ipcw=ipcw,splitMethod=splitMethod,giveToModel=model.args,predictHandlerFun=predictHandlerFun,keep=keep.matrix,verbose=verbose)
    CrossValErr <- kCV$CrossValErr
    if (keep.matrix && B>1)
      CrossValErrMat <- kCV$CrossValErrMat
  }

  # }}}
  # {{{ ----------------------BootstrapCrossValidation----------------------

  if (splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632")){
    if (verbose==TRUE){
      message("Split sample loop (B=",B,")")
    }
    if (missing(testTimes)){
      testTimes <- NULL
    }
    BootCv <- bootstrapCrossValidation(object=object,
                                       data=data,
                                       Y=Y,
                                       status=status,
                                       event=event,
                                       times=times,
                                       cause=cause,
                                       ipcw=ipcw,
                                       ipcw.refit=ipcw.refit,
                                       ipcw.call=ipcw.call,
                                       splitMethod=splitMethod,
                                       multiSplitTest=multiSplitTest,
                                       testIBS=testIBS,
                                       testTimes=testTimes,
                                       confInt=confInt,
                                       confLevel=confLevel,
                                       getFromModel=model.parms,
                                       giveToModel=model.args,
                                       predictHandlerFun=predictHandlerFun,
                                       keepMatrix=keep.matrix,
                                       keepResiduals=keep.residuals,
                                       verbose=verbose,
                                       savePath=savePath,
                                       slaveseed=slaveseed)
    BootstrapCrossValErr <- BootCv$BootstrapCrossValErr
    Residuals <- BootCv$Residuals
    names(BootstrapCrossValErr) <- names(object)
    if (multiSplitTest==TRUE){
      comparisons <- allComparisons(names(object))
      multiSplitTestResults <- list(testIBS=testIBS,B=B,M=M,N=N,testTimes=testTimes)
      multiSplitTestResults$Comparisons <- lapply(1:length(comparisons),function(cc){
        if (length(testTimes)>0){
          allPairwisePvaluesTimes <- do.call("rbind",lapply(BootCv$testedResid,function(b){
            b$pValue[[cc]]}))
          out <- list(pValueTimes=apply(allPairwisePvaluesTimes,2,median))
          if (keep.pvalues==TRUE){
            out$allPairwisePvaluesTimes <- allPairwisePvaluesTimes
          }
        }
        else out <- NULL
        if(length(testIBS)>0){
          allPairwisePvaluesIBS <- sapply(BootCv$testedResid,function(b){
            b$IBSpValue[[cc]]
          })
          out$pValueIBS <- median(allPairwisePvaluesIBS)
        }
        if (keep.pvalues==TRUE){
          out$allPairwisePvaluesIBS <- allPairwisePvaluesIBS}
        out
      })
      names(multiSplitTestResults$Comparisons) <- names(comparisons)
      ## multiSplitTest$splitMethod <- splitMethod
      class(multiSplitTestResults) <- "multiSplitTest"
    }
    ## upperLimits <- lapply(BootCv$testedResid,function(x){x[,1:length(testTimes)]})
    ##     if (testIBS==TRUE){
    ##       wtestIBSpValues <- do.call("cbind",apply(BootCv$testedResid,function(x){x[,length(testTimes)+1]}))
    ##     }
    ## wtestIBSupper <- BootCv$testedResid$wtestIBSupper
    ##   }
    if (keep.matrix==TRUE){
      BootstrapCrossValErrMat <- BootCv$BootstrapCrossValErrMat
      names(BootstrapCrossValErr) <- names(object)
    }
  }

  # }}}
# {{{ Bootstrap .632
  if (splitMethod$internal.name=="Boot632"){
    B632Err <- lapply(1:NF,function(f){
      .368 * AppErr[[f]] + .632 * BootstrapCrossValErr[[f]]
    })
    names(B632Err) <- names(object)
  }
  # }}}    
  # {{{ Bootstrap .632+

  if (splitMethod$internal.name=="Boot632plus"){
    B632plusErr <- lapply(1:NF,function(f){
      Err1 <- pmin(BootstrapCrossValErr[[f]],NoInfErr[[f]])
      overfit <- (Err1 - AppErr[[f]]) / (NoInfErr[[f]] - AppErr[[f]])
      overfit[!(Err1>AppErr[[f]])] <- 0
      w <- .632 / (1 - .368 * overfit)
      B632plusErr <- (1-w) * AppErr[[f]]  + w * Err1
      B632plusErr
      ## w[NoInfErr<=BootstrapCrossValErr] <- 1
      ## B632plus.error <- (1-w) * AppErr  + w * BootstrapCrossValErr
    })
    names(B632plusErr) <- names(object)
  }

  # }}}
  # {{{ prepare output

  out <- switch(splitMethod$internal.name,
                "noPlan"=list("AppErr"=AppErr),
                "Boot632plus"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr,"NoInfErr"=NoInfErr,"Boot632plusErr"=B632plusErr),
                "Boot632"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr,"Boot632Err"=B632Err),
                "BootCv"=list("AppErr"=AppErr,"BootCvErr"=BootstrapCrossValErr),
                "loocv"=list("AppErr"=AppErr,"loocvErr"=CrossValErr),
                "crossval"=list("AppErr"=AppErr,"crossvalErr"=CrossValErr),
                "noinf"=list("AppErr"=AppErr,"NoInfErr"=NoInfErr))
  observed.maxtime <- sapply(out,function(x){
    ## lapply(x,function(y){times[length(y)-sum(is.na(y))-1]})
    lapply(x,function(y){times[length(y)-sum(is.na(y))]})
  })
  minmaxtime <- min(unlist(observed.maxtime))
  if (multiSplitTest==TRUE){
    out <- c(out,list(multiSplitTest=multiSplitTestResults))
  }
  if (keep.residuals==TRUE){
    out <- c(out,list(Residuals=Residuals))
  }
  if (keep.matrix==TRUE && splitMethod$internal.name!="noPlan"){
    if (splitMethod$internal.name %in% c("crossval","loocv")){
      if (B>1)
        out <- c(out,list("CrossValErrMat"=CrossValErrMat))
    }
    else{
      if (splitMethod$internal.name!="noinf")
        out <- c(out,list("BootstrapCrossValErrMat"=BootstrapCrossValErrMat))
    }
  }
  if (!is.na(fillChar))
    out <- lapply(out,function(o){
      o[is.na(o)] <- fillChar
      o
    })
  if (!is.null(model.parms))
    out <- c(out,list("ModelParameters"=BootCv$ModelParameters))
  
  
  if (!keep.index) splitMethod$index <- NULL
  n.risk <- N - sindex(Y,times)
  # }}}
# {{{ put out
  if(keep.models==TRUE)
    outmodels <- object
  else if (keep.models=="Call"){
    outmodels <- lapply(object,function(o){
      cc <- try(o$call,silent=TRUE)
      if(class(cc)=="try-error")
        class(object)
      else
        cc
    })
    names(outmodels) <- names(object)
  }
  else{
    outmodels <- names(object)
    names(outmodels) <- names(object)
  }
  out <- c(out,
           list(call=theCall,
                response=model.response(m),
                time=times,
                ipcw.fit=ipcw$fit,
                n.risk=n.risk,
                models=outmodels,
                maxtime=maxtime,
                observed.maxtime=observed.maxtime,
                minmaxtime=minmaxtime,
                reference=reference,
                start=min(times),
                cens.model=cens.model,
                exact=exact,
                splitMethod=splitMethod))
  ##   if (verbose==TRUE && splitMethod$internal.name %in% c("BootCv","Boot632","Boot632plus","crossval","loocv")) cat("\n")
  class(out) <- "pec"
  out

  # }}}
}


pec.glm <- function(object,...){
  f <- eval(object$call$formula)
  d <- eval(object$call$data)
  pec(list("GLM"=object),formula=f,data=d,...)
}

pec.coxph <- function(object,...){
  f <- eval(object$call$formula)
  d <- eval(object$call$data)
  pec(list("CoxModel"=object),formula=f,data=d,...)
}

pec.cph <- function(object,...){
  f <- eval(object$call$formula)
  d <- eval(object$call$data)
  pec(list("CoxModel"=object),formula=f,data=d,...)
}

pec.survfit <- function(object,...){
  f <- eval(object$call$formula)
  d <- eval(object$call$data)
  pec(list("SurvFit"=object),formula=f,data=d,...)
}

pec.prodlim <- function(object,...){
  f <- eval(object$call$formula)
  d <- eval(object$call$data)
  pec(list("ProductLimit"=object),formula=f,data=d,...)
}

pec.aalen <- function(object,...){
  f <- eval(object$call$formula)
  d <- eval(object$call$data)
  pec(list("AalenModel"=object),formula=f,data=d,...)
}
