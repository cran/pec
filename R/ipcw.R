ipcw <- function(formula,
                 data,
                 model=c("cox","marginal","nonpar","aalen","none"),
                 times,
                 otimes){
  model <- match.arg(model,c("cox","marginal","nonpar","aalen","none"))
  class(model) <- model
  #  print(model)
  UseMethod("ipcw",model)
}

ipcw.none <- function(formula,data,model,times,otimes){
  length.times <- length(times)
  stopifnot(length.times>0)
  wt <- rep(1,length(times))
  wt.obs <- rep(1,length(otimes))
  fit <- NULL
  list(wt=wt,wt.obs=wt.obs,fit=fit)
}


## reverse Kaplan-Meier 
ipcw.marginal <- function(formula,data,model,times,otimes){
  length.times <- length(times)
  stopifnot(length.times>0)
  formula <- update.formula(formula,"~1")
  fit <- prodlim(formula,data=data,reverse=TRUE)
  locotimes <- match(otimes,fit$time,nomatch=NA)
  if (any(is.na(locotimes))) stop("Can not locate all individual observation times" )
  wt.obs <- c(1,fit$surv)[locotimes] ## at (otimes_i-)
  wt <- c(1,fit$surv)[sindex(jump.times=fit$time,eval.times=times) +1] ## at all requested times
  list(wt=wt,wt.obs=wt.obs,fit=fit)
}

## reverse Stone-Beran 
ipcw.nonpar <- function(formula,data,model,times,otimes){
  length.times <- length(times)
  stopifnot(length.times>0)
  fit <- prodlim(formula,data=data,reverse=TRUE,bandwidth="smooth")
  wt <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
  wt.obs <- predictSurvIndividual(fit,lag=1)
  list(wt=wt,wt.obs=wt.obs,fit=fit)
}

#reverse Cox via the survival package
ipcw.coxph <- function(formula,data,model,times,otimes){
  length.times <- length(times)
  stopifnot(length.times>0)
  require(survival)
  survival.survfit.coxph <- getFromNamespace("survfit.coxph",ns="survival")
  survival.summary.survfit <- getFromNamespace("summary.survfit",ns="survival")
  status.name <- all.vars(formula)[2]
  reverse.data <- data
  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  stopifnot(NROW(na.omit(data))>0)
  fit <- coxph(formula,data=reverse.data)
  survfit.object <- survival.survfit.coxph(fit,newdata=data,se.fit=FALSE,conf.int=FALSE)
  wt <- survest(fit,newdata=data,times=times,se.fit=FALSE)$surv
  wt.obs <- survest(fit,times=otimes-min(diff(c(0,unique(otimes))))/2,what='parallel')
  list(wt=wt,wt.obs=wt.obs,fit=fit)
}

#reverse Cox via Harrel's package
ipcw.cox <- function(formula,data,model,times,otimes){
  length.times <- length(times)
  stopifnot(length.times>0)
  require(Design)
  status.name <- all.vars(formula)[2]
  reverse.data <- data
  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  stopifnot(NROW(na.omit(data))>0)
  fit <- cph(formula,data=reverse.data,surv=TRUE,x=TRUE,y=TRUE)
  wt <- survest(fit,newdata=data,times=times,se.fit=FALSE)$surv
  wt.obs <- survest(fit,times=otimes-min(diff(c(0,unique(otimes))))/2,what='parallel')
  list(wt=wt,wt.obs=wt.obs,fit=fit)
}

#reverse Aalen model via the timereg package
ipcw.aalen <- function(formula,data,model,times,otimes){
  length.times <- length(times)
  stopifnot(length.times>0)
  require(timereg)
  require(Design)
  status.name <- all.vars(formula)[2]
  reverse.data <- data
  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
  fit <- do.call(model,list(formula=formula,data=reverse.data,n.sim=0))
  wt <- survest(fit,newdata=data,times=times)
  wt.obs <- diag(survest(fit,newdata=data,times=pmax(0,otimes-min(diff(unique(otimes)))/2)))
  list(wt=wt,wt.obs=wt.obs,fit=fit)
}

#Stone-Beran using the linear predictor of a reverse Cox model
#ipcw.project <- function(formula,data,model,times,otimes){ 
#  length.times <- length(times)
#  stopifnot(length.times>0)
#  require(Design)
#  status.name <- all.vars(formula)[2]
#  reverse.data <- data
#  reverse.data[,status.name] <- 1 * (reverse.data[,status.name]==0)
#  stopifnot(NROW(na.omit(data))>0)
#  fit <- cph(formula,data=reverse.data,surv=TRUE,x=TRUE,y=TRUE)
#  data$lp <- predict(fit,type="lp")
#  reform <- reformulate("lp",formula[[2]])
#  fit <- prodlim(reform,data=data,reverse=TRUE,bandwidth="smooth")
#  wt <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
#  wt.obs <- predictSurvIndividual(fit,lag=1)
#  list(wt=wt,wt.obs=wt.obs,fit=fit)
#}
