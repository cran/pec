oldCindex <- function(formula,data,cens.formula=formula,cens.model="cox",conf.int=NULL,replan="none",M,B,maxtime,verbose=TRUE){

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

  if (missing(data)) stop("Argument data is missing.")

  # censoring model
  # --------------------------------------------------------------------
  
  cens.model <- match.arg(cens.model,c("project","cox","marginal","nonpar","aalen","none"))
  
  # response
  # --------------------------------------------------------------------
  m <- model.frame(formula,data,na.action=na.fail)
  response <- model.response(m)
  if (survp==FALSE && NCOL(response)!=1) stop("Response must be one-dimensional.")
  if (survp==TRUE && NCOL(response)!=2) stop("Survival response must at least consist of two columns: time and status.")

  # predictor
  # --------------------------------------------------------------------
  
  Z <- m[,2,drop=TRUE]
  
  # sort the data 
  # --------------------------------------------------------------------

  if (survp){
    neworder <- order(response[,"time"],-response[,"status"])
    response <- response[neworder,,drop=FALSE]
    Y <- response[,"time"]
    status <- response[,"status"]
    Z <- Z[neworder]
  }
  else{
    cens.model <- "none"
    neworder <- order(response)
    Y <- response[neworder]
    status <- rep(1,length(Y))
    Z <- Z[neworder]
  }
  data <- data[neworder,]
  
  times <- unique(Y)
  N <- length(Y)
  NU <- length(times)

  # find the range of the response 
  # --------------------------------------------------------------------
  
  if (missing(maxtime) || is.null(maxtime))
    maxtime <- times[NU]

  # find weights for right censored data (that are 1 if no censoring) 
  # --------------------------------------------------------------------

  if (missing(cens.formula)) cens.formula <- formula
  if((cens.model %in% c("aalen","cox","nonpar"))){
    if (all(as.numeric(status)==1) || sum(status)==N){
      message("No censored observations: cens.model coerced to \"none\".")
      cens.model <- "none"
    }
    if (length(attr(terms(cens.formula),"factors"))==0){
      if (verbose==TRUE)
        message("No covariates  specified: cens.model coerced to \"marginal\".\n")
      cens.model <- "marginal"}
  }

  #  ipcw <- ipcw(formula=cens.formula,data=data,model=cens.model,times=times,otimes=Y)
  #  wt <- ipcw$wt
  #  wt.obs <- ipcw$wt.obs
  #  stopifnot(length(wt.obs)==N)
  #  wt.dim <- if (cens.model%in% c("marginal","none")) 0 else 1

  #  FIXME: what is the correct weight??? G(T_i) or G(T_i-) or G(T_j-) or ...
  
  weight <- ipcw(formula=cens.formula,data=data,model=cens.model,times=times,otimes=Y)$wt.obs
  weight.lag <- ipcw(formula=cens.formula,data=data,model=cens.model,times=times-min(diff(times))/2,otimes=Y)$wt.obs
  
  # cindex
  # --------------------------------------------------------------------
  
  tindex <- match(Y,times)
  cindex <- .C("cindex",cindex=double(1),as.integer(tindex),as.double(Y),as.integer(status),as.double(weight),as.double(weight.lag),as.double(Z),as.integer(N),as.double(maxtime),as.integer(!is.null(dim(weight))),NAOK=TRUE,package="pec")$cindex

  if (!is.null(conf.int)){
    if (!is.numeric(conf.int) || conf.int> 1  || conf.int<0) {
      warning("Incorrect specification of conf.int (the level for confidence interval). Set to most popular value: 0.95!",immediate.=TRUE,call.=TRUE)
      conf.int <- .95
    }
    boot <- matrix(sapply(1:B,function(b){sort(sample(1:N,replace=TRUE))}),nrow=N,ncol=B)
    boot.cindex <- sapply(1:B,function(b){
      who <- boot[,b,drop=TRUE]
      if (!is.null(dim(weight))) {
        weight.b <- weight[who,]
        weight.lag.b <- weight.lag[who,]}
      else{
        weight.b <- weight
        weight.lag.b <- weight.lag
      }
      cindex.b <- .C("cindex",
                     cindex=double(1),
                     as.integer(tindex[who]),
                     as.double(Y[who]),
                     as.integer(status[who]),
                     as.double(weight.b),
                     as.double(weight.lag.b),
                     as.double(Z[who]),
                     as.integer(N),
                     as.double(maxtime),
                     as.integer(!is.null(dim(weight))),
                     NAOK=TRUE,
                     package="pecDev")$cindex
      cindex.b
    })
    alpha <- (1-conf.int)/2
    conf.bounds <- quantile(boot.cindex - cindex,c(alpha,1-alpha))
    #    print(conf.bounds)
    list(cindex=cindex,lower=cindex-conf.bounds[2],upper=cindex-conf.bounds[1])
  }
  else
    cindex
}

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
#  hindex <- .C("hindex",hindex=double(1),as.double(Y),as.integer(status),as.double(Z),as.integer(N),as.integer(maxtime),NAOK=FALSE,PACKAGE="prodlim")$hindex
#  hindex
#}

