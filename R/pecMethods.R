pec <- function(object,...){
  UseMethod("pec",object=object)
}

pec.list <- function(object,
                     formula,
                     data,
                     times,
                     start,
                     maxtime,
                     exact=TRUE,
                     exactness=100,
                     fillChar=NA,
                     cens.model="cox",
                     replan="none",
                     B,
                     M,
                     model.args=NULL,
                     model.parms=NULL,
                     keep.matrix=FALSE,
                     import=NULL,
                     export=NULL,
                     na.accept=0,
                     verbose=TRUE,
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
    if (match("formula",class(formula),nomatch=0)==0)
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
  
  # find jumptimes in the range of the response 
  # --------------------------------------------------------------------
  
  if (missing(maxtime) || is.null(maxtime))
    maxtime <- unique.Y[NU]
  
  if (missing(times)){
    if (exact==TRUE)
      times <- unique.Y
    else{  
      if (missing(start))
        if (survp==TRUE) start <- 0 ## survival times are positive
        else start <- min(unique.Y)
      times <- seq(start,maxtime,(maxtime - start)/exactness)
    }
  }
  else{
    if (exact==TRUE) 
      times <- sort(c(unique(times),unique.Y))
    else
      times <- sort(unique(times))
  }
  # print(times)
  times <- times[times<=maxtime]
  NT <-  length(times)
  tindex <- sindex(jump.times=unique.Y,eval.times=times)
  
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
  ipcw <- ipcw(formula=formula,data=data,model=cens.model,times=times,otimes=Y)
  wt <- ipcw$wt
  wt.obs <- ipcw$wt.obs
  stopifnot(length(wt.obs)==N)
  wt.dim <- if (cens.model%in% c("marginal","none")) 0 else 1
  
  #  if (NCOL(wt)>1) {stopifnot(length(wt)==(N*NT))}  else{stopifnot(length(wt)==NT)}


  # replan
  # --------------------------------------------------------------------
  replan <- resolveReplan(replan=replan,B=B,N=N,M=M,k=k,import=import,export=export)
  B <- replan$B
  ResampleIndex <- replan$index
  k <- replan$k
  do.resample <- !(is.null(ResampleIndex))

  # checking the models for compatibility with resampling
  # --------------------------------------------------------------------
  if (do.resample){
    cm <- checkModels(object=object,model.args=model.args,model.parms=model.parms,replan=replan$internal.name)
    model.args <- cm$model.args
    model.parms <- cm$model.parms
  }
  
  # computation of prediction error in a loop over the models 
  # --------------------------------------------------------------------
  
  list.out <- lapply(1:NF,function(f){
    
    if (verbose==TRUE) cat("\n",names(object)[f],"\n")
    fit <- object[[f]]
    extract <- model.parms[[f]]
    
    # apparent error (use the same data for fitting and validation)
    # --------------------------------------------------------------------
    
    if (is.null(model.args[[f]])){
      pred <- predictSurvProb(fit,newdata=data,times=times,train.data=data)
    }
    else{
      print(class(fit))
      pred <- do.call("predictSurvProb",c(list(object=fit,newdata=data,times=times,train.data=data),model.args[[f]]))
    }
    
    AppErr <- .C("pec",
                 pec=double(NT),
                 as.double(Y),
                 as.double(status),
                 as.double(times),
                 as.double(pred),
                 as.double(wt),
                 as.double(wt.obs),
                 as.integer(N),
                 as.integer(NT),
                 as.integer(wt.dim),
                 as.integer(NCOL(pred)>1),
                 NAOK=TRUE,
                 PACKAGE="pec")$pec

    # No information error  
    # --------------------------------------------------------------------
    if (replan$internal.name %in% c("boot632plus","noinf")){
      # if (NCOL(pred)==1) NoInfErr <- AppErr
      NoInfErr <- .C("pec_noinf",pec=double(NT),as.double(Y),as.double(status),as.double(times),as.double(pred),as.double(wt),as.double(wt.obs),as.integer(N),as.integer(NT),as.integer(wt.dim),as.integer(NCOL(pred)>1),NAOK=TRUE,PACKAGE="pec")$pec
    }

    # resampling error 
    # --------------------------------------------------------------------
    if (replan$internal.name=="plain"){
      compute.BootErrMat <- lapply(1:B,function(b){
        if (verbose==TRUE) internalTalk(b,B)
        index.b <- ResampleIndex[,b]
        data.b <- data[index.b,,drop=FALSE]
        wt.b <- if (wt.dim==1) wt[index.b,] else wt
        fit.b <- internalReevalFit(object=fit,data=data.b,step=b,silent=na.accept>0,verbose=verbose)
        if (!is.null(extract)) fit.parms <- fit.b[extract]
        else fit.parms <- NULL
        if (is.null(model.args[[f]])){
          pred.b <- predictSurvProb(fit.b,newdata=data.b,times=times,train.data=data.b)
        }
        else{
          pred.b <- do.call("predictSurvProb",c(list(object=fit.b,newdata=data.b,times=times,train.data=data.b),model.args[[f]]))
        }
        innerBootErr <- .C("pec",pec=double(NT),as.double(Y[index.b]),as.double(status[index.b]),as.double(times),as.double(pred.b),as.double(wt.b),as.double(wt.obs[index.b]),as.integer(M),as.integer(NT),as.integer(wt.dim),as.integer(NCOL(pred.b)>1),NAOK=TRUE,PACKAGE="pec")$pec
        list("innerBootErr"=innerBootErr,"fit.parms"=fit.parms)
      })
      BootErrMat <- do.call("cbind",lapply(compute.BootErrMat,function(x)x$innerBootErr))
      if (!is.null(extract))
        fitParms <- lapply(compute.BootErrMat,function(x)x$fit.parms)
      if (na.accept>0)
        BootErr <- apply(BootErrMat,1,function(b) mean(b,na.rm=sum(is.na(b))<na.accept))
      else
        BootErr <- rowMeans(BootErrMat)
    }
    if (length(k)>0){
      # k-fold CrossValidation
      # --------------------------------------------------------------------
      CrossValErrMat <- do.call("cbind",lapply(1:B,function(b){
        if (verbose==TRUE) internalTalk(b,B)
        groups <- ResampleIndex[,b,drop=TRUE]
        ## each subject belongs to exactly one group
        ## the prediction `p[i]' is obtained with the reduced data
        pred.b <- do.call("rbind",lapply(1:k,function(g){
          id <- groups==g
          train.data <- data[!id,,drop=FALSE]
          fit.k <- internalReevalFit(object=fit,data=train.data,step=paste("CV group=",k),silent=na.accept>0,verbose=verbose)
          # fit.k <- with(train.data,eval(f$call))
          val.data <- data[id,,drop=FALSE]
        if (is.null(model.args[[f]])){
            p.group <- predictSurvProb(fit.k,newdata=val.data,times=times,train.data=train.data)
          }
          else{
            p.group <- do.call("predictSurvProb",c(list(object=fit.k,newdata=val.data,times=times,train.data=train.data),model.args[[f]]))
          }
          if(is.null(dim(p.group))) p.group <- do.call("rbind",lapply(1:NROW(val.data),function(x){p.group}))
          p.group
        }))
        pred.b <- pred.b[order(order(groups)),]
        innerCrossValErr <- .C("pec",
                               pec=double(NT),
                               as.double(Y),
                               as.double(status),
                               as.double(times),
                               as.double(pred.b),
                               as.double(wt),
                               as.double(wt.obs),
                               as.integer(N),
                               as.integer(NT),
                               as.integer(wt.dim),
                               as.integer(NCOL(pred.b)>1),
                               NAOK=TRUE,
                               PACKAGE="pec")$pec
        innerCrossValErr
      }))
      if (na.accept>0)
        CrossValErr <- apply(CrossValErrMat,1,function(b) mean(b,na.rm=sum(is.na(b))<na.accept))
      else
        CrossValErr <- rowMeans(CrossValErrMat)
    }
    if (replan$internal.name %in% c("boot632plus","outofbag","boot632")){
      
      # OutOfBagError aka BootstrapCrossValidationError
      # --------------------------------------------------------------------
      compute.OutOfBagErrMat <- lapply(1:B,function(b){
        if (verbose==TRUE) internalTalk(b,B)
        vindex.b <- match(1:N,ResampleIndex[,b],nomatch=0)==0
        val.b <- data[vindex.b,,drop=FALSE]
        train.b <- data[ResampleIndex[,b],,drop=FALSE]
        if (wt.dim==1) wt.b <- wt[vindex.b,] else wt.b <- wt
        fit.b <- internalReevalFit(object=fit,
                                   data=train.b,
                                   step=b,
                                   silent=na.accept>0,
                                   verbose=verbose)
        if (!is.null(extract)) fit.parms <- fit.b[extract]
        else fit.parms <- NULL
        if (is.null(fit.b)){
          failed <- "fit"
          innerOutOfBagErr <- rep(NA,NT)
        }
        else{
          if (is.null(model.args[[f]])){
            try2predict <- try(pred.b <- predictSurvProb(fit.b,newdata=val.b,times=times,train.data=train.b),silent=na.accept>0)
          }
          else{
            try2predict <- try(pred.b <- do.call("predictSurvProb",c(list(object=fit.b,newdata=val.b,times=times,train.data=train.b),model.args[[f]])))
          }
          if (inherits(try2predict,"try-error")==TRUE){
            if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
            failed <- "prediction"
            innerOutOfBagErr <- rep(NA,NT)
          }
          else{
            failed <- NA
            ## if (write.pred==TRUE && file.exists(write.path)){
            ## write.table(pred.b,file=paste(write.path,"pred-oob-sample-",b,"-",names(object)[f],".txt",sep=""))
            ## message(paste("prediction of ",names(object)[f],"for oob-sample",b,"saved"))
            ## }
            innerOutOfBagErr <- .C("pec",pec=double(NT),as.double(Y[vindex.b]),as.double(status[vindex.b]),as.double(times),as.double(pred.b),as.double(wt.b),as.double(wt.obs[vindex.b]),as.integer(sum(vindex.b)),as.integer(NT),as.integer(wt.dim),as.integer(NCOL(pred.b)>1),NAOK=TRUE,PACKAGE="pec")$pec
          }
        }
        list("innerOutOfBagErr"=innerOutOfBagErr,"fit.parms"=fit.parms,"failed"=failed)
      })
      if (verbose==TRUE) cat("\n")
      if (!is.null(extract)) fitParms <- lapply(compute.OutOfBagErrMat,function(x)x$fit.parms)
      failed <- na.omit(sapply(compute.OutOfBagErrMat,function(x)x$failed))
      OutOfBagErrMat <- do.call("cbind",lapply(compute.OutOfBagErrMat,function(x)x$innerOutOfBagErr))
      if (na.accept>0)
        OutOfBagErr <- apply(OutOfBagErrMat,1,function(b) mean(b,na.rm=sum(is.na(b))<na.accept))
      else
        OutOfBagErr <- rowMeans(OutOfBagErrMat)
    }
    # Bootstrap .632
    # --------------------------------------------------------------------
    if (replan$internal.name=="boot632"){
      B632Err <- .368 * AppErr + .632 * OutOfBagErr
    }
    # Bootstrap .632+
    # --------------------------------------------------------------------
    if (replan$internal.name=="boot632plus"){
      Err1 <- pmin(OutOfBagErr,NoInfErr)
      overfit <- (Err1 - AppErr) / (NoInfErr - AppErr)
      overfit[!(Err1>AppErr)] <- 0
      w <- .632 / (1 - .368 * overfit)
      B632plusErr <- (1-w) * AppErr  + w * Err1
      ## w[NoInfErr<=OutOfBagErr] <- 1
      ## B632plus.error <- (1-w) * AppErr  + w * OutOfBagErr
    }
    if (length(k)>0){
      out <- list("PredErr"=CrossValErr,"AppErr"=AppErr)
      if (keep.matrix==TRUE)
        out <- c(out,list("CrossValErrMat"=CrossValErrMat))
    }
    else{
      out <- switch(replan$internal.name,
                    "noPlan"=list("PredErr"=AppErr),
                    "plain"=list("PredErr"=BootErr,"AppErr"=AppErr),
                    "boot632plus"=list("AppErr"=AppErr,"OutOfBagErr"=OutOfBagErr,"NoInfErr"=NoInfErr,"weight"=w,"overfit"=overfit,"PredErr"=B632plusErr),
                    "boot632"=list("AppErr"=AppErr,"OutOfBagErr"=OutOfBagErr,"PredErr"=B632Err),
                    "outofbag"=list("AppErr"=AppErr,"PredErr"=OutOfBagErr),
                    "noinf"=list("AppErr"=AppErr,"PredErr"=NoInfErr))
      if (keep.matrix==TRUE && replan$internal.name!="no"){
        if (replan$internal.name=="plain") out <- c(out,"BootErrMat"=BootErrMat)
        else if (replan$internal.name!="noinf") out <- c(out,list("OutOfBagErrMat"=OutOfBagErrMat))
      }
    }
    if (!is.na(fillChar))
      out <- lapply(out,function(o){
        o[is.na(o)] <- fillChar
        o
      })
    if (!is.null(extract)) out <- c(out,list("fitParms"=fitParms))
    if (na.accept>0) out <- c(out,list("failed"=failed))
    out
  })
  
  names.lout <- names(list.out[[1]])
  out <- lapply(names.lout,function(w){
    e <- lapply(list.out,function(x){x[[w]]})
    names(e) <- names(object)
    e
  })
  names(out) <- names.lout
  
  ## data.frame(matrix(unlist(lapply(list.out,function(x){x[[w]]})),ncol=length(fit.out),nrow=NT,dimnames=list(times,names(object))))
  ##   PredErr <- pec$PredErr
  ##   names(PredErr) <- names(object)
  n.risk <- N - sindex(Y,times)
  
  out <- c(out,
           list(call=match.call(),
                time=times,
                ipcw.fit=ipcw$fit,
                n.risk=n.risk,
                models=object,
                maxtime=maxtime,
                start=min(times),
                cens.model=cens.model,
                exact=exact,
                method=replan))

  ##   if(length(pec)>1) out <- c(out,pec[-match("PredErr",names(pec))])
  if (verbose==TRUE) cat("\n")
  class(out) <- "pec"
  out
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
