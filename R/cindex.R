cindex <- function(object,...){
  UseMethod("cindex",object=object)
}

cindex.list <- function(object,
                        formula,
                        data,
                        eval.times,
                        pred.times,
                        starttime,
                        maxtime,
                        cens.model="marg",
                        splitMethod="noPlan",
                        B,
                        M,
                        model.args=NULL,
                        model.parms=NULL,
                        keepModels=FALSE,
                        keep.weights=FALSE,
                        keep.index=FALSE,
                        keep.matrix=FALSE,
                        ## na.accept=0,
                        verbose=FALSE,
                        ...){

  # models
  # --------------------------------------------------------------------
  NF <- length(object) 
  if (is.null(names(object)))names(object) <- sapply(object,function(o)class(o)[1])
  else{names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])}
  names(object) <- make.names(names(object),unique=TRUE)
  
  # formula
  # --------------------------------------------------------------------

  if (missing(formula)){
    formula <- eval(object[[1]]$call$formula)
    if (class(formula)!="formula")
      stop("Argument formula is missing.")
    else if (verbose)
      warning("Argument formula is missing. I use the formula from the call to the first model instead.")
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
  
  # data
  # --------------------------------------------------------------------

  ##   print(str(eval(object[[1]]$call$data)))
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (match("data.frame",class(data),nomatch=0)==0)
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I use the data from the call to the first model instead.")
  }

  # censoring model
  # --------------------------------------------------------------------
  
  cens.model <- match.arg(cens.model,c("cox","marginal","nonpar","aalen","none"))
  
  # response
  # --------------------------------------------------------------------
  m <- model.frame(formula,data,na.action=na.fail)
  response <- model.response(m)
  if (survp==FALSE && NCOL(response)!=1) stop("Response must be one-dimensional.")
  if (survp==TRUE && NCOL(response)!=2) stop("Survival response must at least consist of two columns: time and status.")
  
  # sort the data 
  # --------------------------------------------------------------------
  if (survp){
    neworder <- order(response[,"time"],-response[,"status"])
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
  data <- data[neworder,]
  
  unique.Y <- unique(Y)
  N <- length(Y)
  NU <- length(unique.Y)
  
  # define the prediction time(s) and the evaluation time(s)
  # --------------------------------------------------------------------
  if (missing(starttime) || is.null(starttime)) starttime <- 0
  if (missing(maxtime) || is.null(maxtime)) maxtime <- unique.Y[NU]+1
  if (missing(eval.times)) {
    eval.times <- max(Y[status==1])
  }
  else{
    eval.times <- eval.times[eval.times<=maxtime]
  }
  
  if (missing(pred.times))
    pred.times <- median(unique.Y)
  
  ## FIXME: if the model changes the risk order over time, then we need to care about pred.times
  
  NT <-  length(eval.times)
  tindex <- match(Y,unique.Y)
  
  # find weights for right censored data (that are 1 if no censoring) 
  # --------------------------------------------------------------------

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
  #  FIXME: what are the correct weights??? G(T_i) or G(T_i-) or G(T_j-) or ...
  #  FIXED: the correct weights are G(T_i|X_?) and G(T_i-|X_?)
  weight.i <- ipcw(formula=formula,
                   data=data,
                   method=cens.model,
                   times=unique.Y,
                   subjectTimes=Y,
                   subjectTimesLag=1,
                   what="IPCW.subjectTimes")$IPCW.subjectTimes
  weight.j <- ipcw(formula=formula,
                   data=data,
                   method=cens.model,
                   times=unique.Y,
                   subjectTimes=Y,
                   subjectTimesLag=0,
                   what="IPCW.times")$IPCW.times

  #  print(dim(weight.i))
  #  print(dim(weight.j))
  #  stop()
  # method
  # --------------------------------------------------------------------
  
  splitMethod <- resolvesplitMethod(splitMethod=splitMethod,B=B,N=N,M=M)
  
  B <- splitMethod$B
  ResampleIndex <- splitMethod$index
  if (!keep.index) splitMethod$index <- NULL
  k <- splitMethod$k
  do.resample <- !(is.null(ResampleIndex))
  
  # checking the models for compatibility with resampling
  # --------------------------------------------------------------------
  if (do.resample){
    cm <- checkModels(object=object,
                      model.args=model.args,
                      model.parms=model.parms,
                      splitMethod=splitMethod$internal.name)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }
  
  # computation of the cindex in a loop over the models 
  # --------------------------------------------------------------------

  list.out <- lapply(1:NF,function(f){
    if (verbose==TRUE) cat("\n",names(object)[f],"\n")
    fit <- object[[f]]
    extract <- model.parms[[f]]
    
    # apparent error (use the same data for fitting and validation)
    # --------------------------------------------------------------------
    #    otimes <- Y
    pred <- do.call("predictSurvProb",
                    c(list(object=fit,
                           newdata=data,
                           times=pred.times,
                           train.data=data),
                      model.args[[f]]))
    
    # print(cbind(Y,status,pred))
    # stop()
    if (length(pred.times)==1 && length(pred.times)<length(eval.times))
      pred <- rep(pred,length(eval.times))
    AppCindexResult <- .C("cindex",cindex=double(NT),conc=double(NT),pairs=double(NT),as.integer(tindex),as.double(Y),as.integer(status),as.double(eval.times),as.double(weight.i),as.double(weight.j),as.double(pred),as.integer(N),as.integer(NT),as.double(starttime),as.double(maxtime),as.integer(!is.null(dim(weight.j))),NAOK=TRUE,package="pec")
    AppCindex <- AppCindexResult$cindex
    AppPairs <- AppCindexResult$pairs
    AppConcordant <- AppCindexResult$conc
    if (splitMethod$internal.name %in% c("Boot632plus","BootCv","Boot632")){
      # BootCvError or BootstrapCrossValidationCindexor
      # --------------------------------------------------------------------
      compute.BootCvCindexList <- lapply(1:B,function(b){
        if (verbose==TRUE) internalTalk(b,B)
        vindex.b <- match(1:N,unique(ResampleIndex[,b]),nomatch=0)==0
        val.b <- data[vindex.b,,drop=FALSE]
        train.b <- data[ResampleIndex[,b],,drop=FALSE]
        fit.b <- internalReevalFit(object=fit,data=train.b,step=b,silent=FALSE,verbose=verbose)
        if (!is.null(extract)) fit.parms <- fit.b[extract]
        else fit.parms <- NULL
        if (is.null(fit.b)){
          failed <- "fit"
          innerBootCvCindex <- NA
          innerBootCvPairs <- NA
          innerBootCvConcordant <- NA
        }
        else{
          try2predict <- try(pred.b <- do.call("predictSurvProb",c(list(object=fit.b,newdata=val.b,times=pred.times,train.data=train.b),model.args[[f]])))
          if (inherits(try2predict,"try-error")==TRUE){
            if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
            failed <- "prediction"
            innerBootCvCindex <- NA
            innerBootCvPairs <- NA
            innerBootCvConcordant <- NA
          }
          else{
            failed <- NA
            Y.b <- Y[vindex.b]
            tindex.b <- match(Y.b,unique(Y.b))
            innerBootCvLoop <- .C("cindex",cindex=double(NT),conc=double(NT),pairs=double(NT),as.integer(tindex.b),as.double(Y.b),as.integer(status[vindex.b]),as.double(eval.times),as.double(weight.i[vindex.b]),as.double(weight.j[vindex.b]),as.double(pred.b),as.integer(sum(vindex.b)),as.integer(NT),as.double(starttime),as.double(maxtime),as.integer(!is.null(dim(weight.j))),NAOK=TRUE,package="pec")
            innerBootCvCindex <- innerBootCvLoop$cindex
            innerBootCvPairs <- innerBootCvLoop$pairs
            innerBootCvConcordant <- innerBootCvLoop$conc
          }
        }
        list("innerBootCvConcordant"=innerBootCvConcordant,
             "innerBootCvPairs"=innerBootCvPairs,
             "innerBootCvCindex"=innerBootCvCindex,
             "fit.parms"=fit.parms,
             "failed"=failed)
      })
      if (verbose==TRUE) cat("\n")
      if (!is.null(extract)) fitParms <- lapply(compute.BootCvCindexList,function(x)x$fit.parms)
      failed <- na.omit(sapply(compute.BootCvCindexList,function(x)x$failed))
      BootCvCindexList <- do.call("cbind",lapply(compute.BootCvCindexList,function(x)x$innerBootCvCindex))
      BootCvConcordantList <- do.call("cbind",lapply(compute.BootCvCindexList,function(x)x$innerBootCvConcordant))
      BootCvPairsList <- do.call("cbind",lapply(compute.BootCvCindexList,function(x)x$innerBootCvPairs))
      ##       if (na.accept>0)
      ##         BootCvCindex <- apply(BootCvCindexList,1,function(b) mean(b,na.rm=sum(is.na(b))<na.accept))
      ##       else
      BootCvCindex <- rowMeans(BootCvCindexList)
      ##       if (na.accept>0)
      ##         BootCvConcordant <- apply(BootCvConcordantList,1,function(b) mean(b,na.rm=sum(is.na(b))<na.accept))
      ##       else
      BootCvConcordant <- rowMeans(BootCvConcordantList)
      ##       if (na.accept>0)
      ##         BootCvPairs <- apply(BootCvPairsList,1,function(b) mean(b,na.rm=sum(is.na(b))<na.accept))
      ##       else
      BootCvPairs <- rowMeans(BootCvPairsList)
    }
    # Bootstrap .632
    # --------------------------------------------------------------------
    if (splitMethod$internal.name=="Boot632"){
      B632Cindex <- .368 * AppCindex + .632 * BootCvCindex
    }
    out <- switch(splitMethod$internal.name,
                  "noPlan"=list("PredCindex"=AppCindex,
                    "Pairs"=AppPairs,
                    "Concordant"=AppConcordant),
                  "Boot632"=list("AppCindex"=AppCindex,
                    "Pairs"=AppPairs,
                    "Concordant"=AppConcordant,
                    "BootCvCindex"=BootCvCindex,
                    "PredCindex"=B632Cindex),
                  "BootCv"=list("AppCindex"=AppCindex,
                    "Pairs"=AppPairs,
                    "Concordant"=AppConcordant,
                    "PredCindex"=BootCvCindex,
                    "BootCvConcordant"=BootCvConcordant,
                    "BootCvPairs"=BootCvPairs
                    ))
    #    if (keep.pred==TRUE && splitMethod$internal.name!="no"){
    #      out <- c(out,list("BootCvCindexList"=BootCvCindexList))
    #    }
    if (keep.matrix==TRUE && splitMethod$internal.name!="no"){
      out <- c(out,list("BootCvCindexList"=BootCvCindexList))
    }
    if (!is.null(extract)) out <- c(out,list("fitParms"=fitParms))
    ## if (na.accept>0) out <- c(out,list("failed"=failed))
    out
  })
  names.lout <- names(list.out[[1]])
  out <- lapply(names.lout,function(w){
    e <- lapply(list.out,function(x){x[[w]]})
    names(e) <- names(object)
    e
  })
  names(out) <- names.lout
  if(keepModels==TRUE)
    outmodels <- object
  else if (keepModels=="Call"){
    outmodels <- lapply(object,function(o)o$call)
    names(outmodels) <- names(object)
  }
  else{
    outmodels <- names(object)
    names(outmodels) <- names(object)
  }
  if(keep.weights)
    www <- list(weight.i,weight.j)
  else
    {www <- NULL
   }
  out <- c(out,
           list(call=match.call(),
                time=eval.times,
                pred.time=pred.times,
                models=outmodels,
                splitMethod=splitMethod,
                weights=www,
                cens.model=cens.model,
                starttime=starttime,
                maxtime=maxtime))
  if (verbose==TRUE) cat("\n")
  class(out) <- "Cindex"
  out
}
    
  
