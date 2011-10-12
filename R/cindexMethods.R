cindex <- function(object,...){
  UseMethod("cindex",object=object)
}
# {{{ header cindex.list

cindex.list <- function(object,
                        formula,
                        data,
                        eval.times,
                        pred.times,
                        cause,
                        cens.model="marginal",
                        ipcw.refit=FALSE,
                        tiedPredictionsIn=TRUE,
                        tiedOutcomeIn=TRUE,
                        tiedMatchIn=TRUE,
                        splitMethod="noPlan",
                        B,
                        M,
                        model.args=NULL,
                        model.parms=NULL,
                        keep.models="Call",
                        keep.residuals=FALSE,
                        keep.pvalues=FALSE,
                        keep.weights=FALSE,
                        keep.index=FALSE,
                        keep.matrix=FALSE,
                        multiSplitTest=FALSE,
                        testTimes,
                        confInt=FALSE,
                        confLevel=0.95,
                        verbose=TRUE,
                        savePath=NULL,
                        ...){

  # }}}
  # {{{ checking integrity some arguments
  theCall=match.call()
  if (match("replan",names(theCall),nomatch=FALSE))
    stop("Argument name 'replan' has been replaced by 'splitMethod'.")
  
  if (keep.residuals && missing(testTimes))
    stop("To keep.residuals please specify testTimes.")
  if (missing(splitMethod) && multiSplitTest==TRUE){
    stop("Need data splitting to compute van de Wiel's test")
  }
  if (missing(M) && multiSplitTest) M <- NA
  stopifnot(as.numeric(tiedPredictionsIn) %in% c(0,1))
  stopifnot(as.numeric(tiedOutcomeIn) %in% c(0,1))
  stopifnot(as.numeric(tiedMatchIn) %in% c(0,1))
  # }}}
  # {{{  formula
  if (missing(formula)){
    formula <- eval(object[[1]]$call$formula)
    if (class(formula)!="formula")
      stop("Argument formula is missing.")
    else if (verbose)
      warning("Formula missing. Using formula from first model")
  }
  formula.names <- try(all.names(formula),silent=TRUE)
  if (!(formula.names[1]=="~")
      ||
      (match("$",formula.names,nomatch=0)+match("[",formula.names,nomatch=0)>0)){
    stop("Invalid specification of formula. Perhaps forgotten right hand side?\nNote that any subsetting, ie data$var or data[,\"var\"], is invalid for this function.")
  }
  else{
    if (!(formula.names[2] %in% c("Surv","Hist")))
      survp <- FALSE
    else
      survp <- TRUE
  }
  # }}}
  # {{{  data
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (match("data.frame",class(data),nomatch=0)==0)
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I use the data from the call to the first model instead.")
  }

# }}}
  # {{{  censoring model
  
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
  if (predictHandlerFun=="predictEventProb")
    stop("Not yet defined: cindex for competing risks")
  # }}}
  # {{{ prediction models
  NF <- length(object) 
  if (is.null(names(object)))names(object) <- sapply(object,function(o)class(o)[1])
  else{names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])}
  names(object) <- make.names(names(object),unique=TRUE)
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
  else{
    event <- NULL
  }
  data <- data[neworder,]
  unique.Y <- unique(Y)
  N <- length(Y)
  NU <- length(unique.Y)

  # }}}
  # {{{  splitMethod
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
  # {{{  define the prediction time(s) and the evaluation time(s)
  maxtime <- unique.Y[NU] 
  if (missing(eval.times) || is.infinite(eval.times)){
    ## eval.times <- max(Y) ## maybe less efficient 
    eval.times <- max(Y[status==1])
  }
  else{
    tooLate <- sum(eval.times>=maxtime)
    if (tooLate>0){
      if (verbose)
        warning(tooLate," eval.times beyond the maximal evaluation time: ",ifelse(maxtime>1,round(maxtime,1),round(maxtime,3)))
      eval.times <- c(eval.times[eval.times<maxtime],maxtime)
    }
  }
  if (missing(pred.times))
    pred.times <- median(unique.Y)
  ## FIXME: if the model changes the risk order over time, then we need to care about pred.times
  NT <-  length(eval.times)
  tindex <- match(Y,unique.Y)
  # }}}
  # {{{  IPCW (all equal to 1 without censoring)
  if((cens.model %in% c("aalen","cox","nonpar"))){
    if (all(as.numeric(status)==1) || sum(status)==N){
      message("No censored observations: cens.model coerced to \"none\".")
      cens.model <- "none"
    }
    if (length(attr(terms(formula),"factors"))==0){
      if (verbose==TRUE)
        message("No covariates  specified: cens.model coerced to \"marginal\".\n")
      cens.model <- "marginal"}
  }
  #  weights for T_i<=T_j
  #  FIXME: what are the correct weights??? 
  #  FIXED: the correct weights are G(T_i|X_j) and G(T_i-|X_i)
  if (ipcw.refit==TRUE && splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632"))
    ipcw.call <- list(weight.i=list(formula=formula,data=data,method=cens.model,times=unique.Y,subjectTimes=Y,subjectTimesLag=1,what="IPCW.subjectTimes"),
                      weight.j=list(formula=formula,data=data,method=cens.model,times=unique.Y,subjectTimes=Y,subjectTimesLag=0,what="IPCW.times"))
  else
    ipcw.call <- NULL
  weight.i <- ipcw(formula=formula,data=data,method=cens.model,times=unique.Y,subjectTimes=Y,subjectTimesLag=1,what="IPCW.subjectTimes")$IPCW.subjectTimes
  weight.j <- ipcw(formula=formula,data=data,method=cens.model,times=unique.Y,subjectTimes=Y,subjectTimesLag=0,what="IPCW.times")$IPCW.times
  weights <- list(weight.i=weight.i,weight.j=weight.j)
  # }}}
  # {{{  checking the models for compatibility with resampling
  if (do.resample){
    cm <- checkModels(object=object,model.args=model.args,model.parms=model.parms,splitMethod=splitMethod$internal.name)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }
  # }}}
  # {{{ -------------------Apparent or test sample cindex----------------------
  
  AppCindexList <- lapply(1:NF,function(f){
    fit <- object[[f]]
    extraArgs <- model.args[[f]]
    pred <- do.call(predictHandlerFun,c(list(object=fit,newdata=data,times=pred.times,train.data=data),extraArgs))
    if (length(pred.times)==1 && length(pred.times)<length(eval.times))
      pred <- rep(pred,length(eval.times))
    if (predictHandlerFun=="predictEventProb"){
      stop("Not yet defined: cindex for competing risks")
    }
    else{
      AppCindexResult <- .C("cindex",
                            cindex=double(NT),
                            conc=double(NT),
                            pairs=double(NT),
                            as.integer(tindex),
                            as.double(Y),
                            as.integer(status),
                            as.double(eval.times),
                            as.double(weight.i),
                            as.double(weight.j),
                            as.double(pred),
                            as.integer(N),
                            as.integer(NT),
                            as.integer(tiedPredictionsIn),
                            as.integer(tiedOutcomeIn),
                            as.integer(tiedMatchIn),
                            as.integer(!is.null(dim(weight.j))),
                            NAOK=TRUE,
                            package="pec")
      AppCindex <- AppCindexResult$cindex
      AppPairs <- AppCindexResult$pairs
      AppConcordant <- AppCindexResult$conc
      list(AppCindex=AppCindex,AppPairs=AppPairs,AppConcordant=AppConcordant)
    }
  })
  AppCindex <- lapply(AppCindexList,function(x){
    x$AppCindex
  })
  AppPairs <- lapply(AppCindexList,function(x){
    2*x$AppPairs
  })
  AppConcordant <- lapply(AppCindexList,function(x){
    2*x$AppConcordant
  })
  names(AppCindex) <- names(object)
  names(AppPairs) <- names(object)
  names(AppConcordant) <- names(object)

  # }}}
  # {{{ ----------------------BootstrapCrossValidation----------------------
  if (splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632")){
    if (missing(testTimes)){
      testTimes <- NULL
    }
    BootCv <- CindexBootstrapCrossValidation(object=object,
                                             data=data,
                                             Y=Y,
                                             status=status,
                                             event=event,
                                             eval.times=eval.times,
                                             pred.times=pred.times,
                                             weights=weights,
                                             ipcw.refit=ipcw.refit,
                                             ipcw.call=ipcw.call,
                                             splitMethod=splitMethod,
                                             multiSplitTest=multiSplitTest,
                                             testTimes=testTimes,
                                             confInt=confInt,
                                             confLevel=confLevel,
                                             getFromModel=model.parms,
                                             giveToModel=model.args,
                                             predictHandlerFun=predictHandlerFun,
                                             tiedPredictionsIn=tiedPredictionsIn,
                                             tiedOutcomeIn=tiedOutcomeIn,
                                             tiedMatchIn=tiedMatchIn,
                                             keepMatrix=keep.matrix,
                                             keepResiduals=keep.residuals,
                                             verbose=verbose,
                                             savePath=savePath)
    BootstrapCrossValCindex <- BootCv$BootstrapCrossValCindex
    Residuals <- BootCv$Residuals
    names(BootstrapCrossValCindex) <- names(object)
    if (multiSplitTest==TRUE){
      comparisons <- allComparisons(names(object))
      multiSplitTest <- list(B=B,M=M,N=N,testTimes=testTimes)
      multiSplitTest$Comparisons <- lapply(1:length(comparisons),function(cc){
        if (length(testTimes)>0){
          allPairwisePvaluesTimes <- do.call("rbind",lapply(BootCv$testedResid,function(b){
            b$pValue[[cc]]}))
          out <- list(pValueTimes=apply(allPairwisePvaluesTimes,2,median))
          if (keep.pvalues==TRUE){
            out$allPairwisePvaluesTimes <- allPairwisePvaluesTimes
          }
        }
        else out <- NULL
        ## if (keep.pvalues==TRUE){
        ## out$allPairwisePvaluesIBS <- allPairwisePvaluesIBS}
        out
      })
      names(multiSplitTest$Comparisons) <- names(comparisons)
      class(multiSplitTest) <- "multiSplitTest"
    }
    if (keep.matrix==TRUE){
      BootstrapCrossValCindexMat <- BootCv$BootstrapCrossValCindexMat
      names(BootstrapCrossValCindex) <- names(object)
    }
  }

  # }}}
  # {{{ Bootstrap .632
  if (splitMethod$internal.name=="Boot632"){
    B632Cindex <- lapply(1:NF,function(f){
      .368 * AppCindex[[f]] + .632 * BootstrapCrossValCindex[[f]]
    })
    names(B632Cindex) <- names(object)
  }
  # }}}    
  # {{{ prepare output
  out <- switch(splitMethod$internal.name,
                  "noPlan"=list("AppCindex"=AppCindex,
                    "Pairs"=AppPairs,
                    "Concordant"=AppConcordant),
                "Boot632"=list("AppCindex"=AppCindex,
                  ## "Pairs"=AppPairs,
                  ## "Concordant"=AppConcordant,
                  "BootCvCindex"= BootstrapCrossValCindex,
                  "B632Cindex"=B632Cindex),
                "BootCv"=list("AppCindex"=AppCindex,
                  ## "Pairs"=AppPairs,
                  ## "Concordant"=AppConcordant,
                  "BootCvCindex"=BootstrapCrossValCindex
                  ## "BootCvConcordant"=BootCvConcordant,
                  ## "BootCvPairs"=BootCvPairs
                  ))
  observed.maxtime <- sapply(out,function(x){
    lapply(x,function(y){
      eval.times[length(y)-sum(is.na(y))]
    })
  })
  minmaxtime <- min(unlist(observed.maxtime))
  
  
  if (multiSplitTest==TRUE){
    out <- c(out,list(multiSplitTest=multiSplitTest))
  }
  if (keep.residuals==TRUE){
    out <- c(out,list(Residuals=Residuals))
  }
  if (keep.matrix==TRUE && splitMethod$internal.name!="noPlan"){
    if (splitMethod$internal.name %in% c("crossval","loocv")){
      if (B>1)
        out <- c(out,list("BootstrapCrossValCindexMat"=BootstrapCrossValCindexMat))
    }
    else{
      if (splitMethod$internal.name!="noinf")
        out <- c(out,list("BootstrapCrossValCindexMat"=BootstrapCrossValCindexMat))
    }
  }
  if (!is.null(model.parms)) out <- c(out,list("ModelParameters"=BootCv$ModelParameters))
  if (!keep.index) splitMethod$index <- NULL
  if(keep.models==TRUE)
    outmodels <- object
  else if (keep.models=="Call"){
    outmodels <- lapply(object,function(o){
      cc <- try(o$call,silent=TRUE)
      if(class(cc)=="try-error")
        class(o)
      else
        cc
    })
    names(outmodels) <- names(object)
  }
  else{
    outmodels <- names(object)
    names(outmodels) <- names(object)
  }
  
  out <- c(out,list(call=theCall,
                    time=eval.times,
                    pred.time=pred.times,
                    response=model.response(m),
                    models=outmodels,
                    splitMethod=splitMethod,
                    weights=weights,
                    cens.model=cens.model,
                    minmaxtime=minmaxtime,
                    maxtime=maxtime))
  if (verbose==TRUE && do.resample==TRUE) cat("\n")
  # }}}
  class(out) <- "Cindex"
  out
}
    
  
