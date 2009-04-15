cindex <- function(object,...){
  UseMethod("cindex",object=object)
}

cindex.list <- function(object,
                        formula,
                        data,
                        times,
                        maxtime,
                        cens.model="marg",
                        replan="noPlan",
                        B,
                        M,
                        model.args=NULL,
                        model.parms=NULL,
                        keep.matrix=FALSE,
                        na.accept=0,
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
  
  # find the range of the response 
  # --------------------------------------------------------------------
  
  if (missing(maxtime) || is.null(maxtime))
    maxtime <- unique.Y[NU]

  ## FIXME: does the model change the risk order over time?
  if (missing(times)){
    sometime <- round(median(1:NU))
    times <- unique.Y[c(sometime,sometime+1)]
  }    
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

  #  FIXME: what are the correct weights??? G(T_i) or G(T_i-) or G(T_j-) or ...
  
  weight <- ipcw(formula=formula,data=data,model=cens.model,times=unique.Y,otimes=Y)$wt.obs
  weight.lag <- ipcw(formula=formula,data=data,model=cens.model,times=unique.Y-min(diff(unique.Y))/2,otimes=Y)$wt.obs
  
  # replan
  # --------------------------------------------------------------------
  replan <- resolveReplan(replan=replan,B=B,N=N,M=M,k=k,import=NULL,export=NULL)

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
  
  # computation of the cindex in a loop over the models 
  # --------------------------------------------------------------------

  list.out <- lapply(1:NF,function(f){
    if (verbose==TRUE) cat("\n",names(object)[f],"\n")
    fit <- object[[f]]
    ##     print(fit)
    extract <- model.parms[[f]]
    
    # apparent error (use the same data for fitting and validation)
    # --------------------------------------------------------------------
    
    if (is.null(model.args[[f]])){
      pred <- predictSurvOrder(fit,newdata=data,times=times,train.data=data)
    }
    else{
##       print(class(fit))
      pred <- do.call("predictSurvOrder",c(list(object=fit,newdata=data,times=times,train.data=data),model.args[[f]]))
    }
    
    AppCindex <- .C("cindex",cindex=double(1),as.integer(tindex),as.double(Y),as.integer(status),as.double(weight),as.double(weight.lag),as.double(pred),as.integer(N),as.double(maxtime),as.integer(!is.null(dim(weight))),NAOK=TRUE,package="pec")$cindex

    #    print(AppCindex)
    
    if (replan$internal.name %in% c("boot632plus","outofbag","boot632")){
      
      # OutOfBagError or BootstrapCrossValidationCindexor
      # --------------------------------------------------------------------
      compute.OutOfBagCindexList <- lapply(1:B,function(b){
        if (verbose==TRUE) internalTalk(b,B)
        vindex.b <- match(1:N,ResampleIndex[,b],nomatch=0)==0
        val.b <- data[vindex.b,,drop=FALSE]
        train.b <- data[ResampleIndex[,b],,drop=FALSE]
        fit.b <- internalReevalFit(object=fit,data=train.b,step=b,silent=na.accept>0,verbose=verbose)
        if (!is.null(extract)) fit.parms <- fit.b[extract]
        else fit.parms <- NULL
        if (is.null(fit.b)){
          failed <- "fit"
          innerOutOfBagCindex <- NA
        }
        else{
          if (is.null(model.args[[f]])){
            try2predict <- try(pred.b <- predictSurvOrder(fit.b,newdata=val.b,times=times,train.data=train.b),silent=na.accept>0)
          }
          else{
            try2predict <- try(pred.b <- do.call("predictSurvOrder",c(list(object=fit.b,newdata=val.b,times=times,train.data=train.b),model.args[[f]])))
          }
          if (inherits(try2predict,"try-error")==TRUE){
            if (verbose==TRUE) warning(paste("During bootstrapping: prediction for model ",class(fit.b)," failed in step ",b),immediate.=TRUE)
            failed <- "prediction"
            innerOutOfBagCindex <- NA
          }
          else{
            failed <- NA
            innerOutOfBagCindex <- .C("cindex",cindex=double(1),as.integer(tindex),as.double(Y[vindex.b]),as.integer(status[vindex.b]),as.double(weight[vindex.b]),as.double(weight.lag[vindex.b]),as.double(pred.b),as.integer(N),as.double(maxtime),as.integer(!is.null(dim(weight))),NAOK=TRUE,package="pec")$cindex
            # print(innerOutOfBagCindex)
          }
        }
        list("innerOutOfBagCindex"=innerOutOfBagCindex,"fit.parms"=fit.parms,"failed"=failed)
      })
      if (verbose==TRUE) cat("\n")
      if (!is.null(extract)) fitParms <- lapply(compute.OutOfBagCindexList,function(x)x$fit.parms)
      failed <- na.omit(sapply(compute.OutOfBagCindexList,function(x)x$failed))
      OutOfBagCindexList <- do.call("cbind",lapply(compute.OutOfBagCindexList,function(x)x$innerOutOfBagCindex))
      if (na.accept>0)
        OutOfBagCindex <- apply(OutOfBagCindexList,1,function(b) mean(b,na.rm=sum(is.na(b))<na.accept))
      else
        OutOfBagCindex <- rowMeans(OutOfBagCindexList)
    }
    # Bootstrap .632
    # --------------------------------------------------------------------
    if (replan$internal.name=="boot632"){
      B632Cindex <- .368 * AppCindex + .632 * OutOfBagCindex
    }
    # Bootstrap .632+
    # --------------------------------------------------------------------
    if (replan$internal.name=="boot632plus"){
      stop("No .632+ estimate of the c-index available.")
      #      Cindex1 <- pmin(OutOfBagCindex,NoInfCindex)
      #      overfit <- (Cindex1 - AppCindex) / (NoInfCindex - AppCindex)
      #      overfit[!(Cindex1>AppCindex)] <- 0
      #      w <- .632 / (1 - .368 * overfit)
      #      B632plusCindex <- (1-w) * AppCindex  + w * Cindex1
    }
    out <- switch(replan$internal.name,
                  "noPlan"=list("PredCindex"=AppCindex),
                  # "plain"=list("PredCindex"=BootCindex,"AppCindex"=AppCindex),
                  # "boot632plus"=list("AppCindex"=AppCindex,"OutOfBagCindex"=OutOfBagCindex,"NoInfCindex"=NoInfCindex,"weight"=w,"overfit"=overfit,"PredCindex"=B632plusCindex),
                  "boot632"=list("AppCindex"=AppCindex,"OutOfBagCindex"=OutOfBagCindex,"PredCindex"=B632Cindex),
                  "outofbag"=list("AppCindex"=AppCindex,"PredCindex"=OutOfBagCindex),
                  # "noinf"=list("AppCindex"=AppCindex,"PredCindex"=NoInfCindex)
                  )
    if (keep.matrix==TRUE && replan$internal.name!="no"){
      #      if (replan$internal.name=="plain")
      #        out <- c(out,"BootCindexList"=BootCindexList)
      #      else
      #      if (replan$internal.name!="noinf")
      out <- c(out,list("OutOfBagCindexList"=OutOfBagCindexList))
    }
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
  out <- c(out,
           list(call=match.call(),
                models=object,
                method=replan))
  if (verbose==TRUE) cat("\n")
  class(out) <- "Cindex"
  out
}
    
  
  #  if (!is.null(conf.int)){
  #    if (!is.numeric(conf.int) || conf.int> 1  || conf.int<0) {
  #      warning("Incorrect specification of conf.int (the level for confidence interval). Set to most popular value: 0.95!",immediate.=TRUE,call.=TRUE)
  #      conf.int <- .95
  #    }
  #    boot <- matrix(sapply(1:B,function(b){sort(sample(1:N,replace=TRUE))}),nrow=N,ncol=B)
  #    boot.cindex <- sapply(1:B,function(b){
  #      who <- boot[,b,drop=TRUE]
  #      if (!is.null(dim(weight))) {
  #        weight.b <- weight[who,]
  #        weight.lag.b <- weight.lag[who,]}
  #      else{
  #        weight.b <- weight
  #        weight.lag.b <- weight.lag
  #      }
  #      cindex.b <- .C("cindex",
  #                     cindex=double(1),
  #                     as.integer(tindex[who]),
  #                     as.double(Y[who]),
  #                     as.integer(status[who]),
  #                     as.double(weight.b),
  #                     as.double(weight.lag.b),
  #                     as.double(Z[who]),
  #                     as.integer(N),
  #                     as.double(maxtime),
  #                     as.integer(!is.null(dim(weight))),
  #                     NAOK=TRUE,
  #                     package="pec")$cindex
  #      cindex.b
  #    })
  #    alpha <- (1-conf.int)/2
  #    conf.bounds <- quantile(boot.cindex - cindex,c(alpha,1-alpha))
  #    #    print(conf.bounds)
  #    list(cindex=cindex,lower=cindex-conf.bounds[2],upper=cindex-conf.bounds[1])
  #  }
  #  else
#  cindex
#}

#harrell <- function(formula,
#                    data,
#                    maxtime=0){
#  m <- surv.model.frame(formula,data)
#  survobject <- model.response(m)
#  data <- data[attr(m,"survorder"),]
#  Y <- survobject[,1] 
#  status <- survobject[,2]
#  times <- unique(Y)
#  N <- length(Y)
#  Z <- m[,2,drop=TRUE]
#  hindex <- .C("hindex",hindex=double(1),as.double(Y),as.integer(status),as.double(Z),as.integer(N),as.integer(maxtime),NAOK=FALSE,PACKAGE="pec")$hindex
#  hindex
#}

