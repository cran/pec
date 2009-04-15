predictSurvOrder <- function(object,newdata,times,...){
  UseMethod("predictSurvOrder",object)
}

predictSurvOrder.default <- function(object,newdata,times,...){
  stop("No method for evaluating predicted probabilities from objects in class: ",class(object),call.=FALSE)
}

predictSurvOrder.numeric <- function(object,newdata,times,...){
  if (!is.numeric(object))
    p <- as.numeric(object)
  else
    p <- object
  if (NROW(p) != NROW(newdata))
    stop("No method for evaluating predicted probabilities from objects in class: ",class(object),call.=FALSE)
  else
  p
}

predictSurvOrder.vector <- function(object,newdata,times,...){
  if (!is.numeric(object))
    p <- as.numeric(object)
  else
    p <- object
  if (NROW(p) != NROW(newdata))
    stop("No method for evaluating predicted probabilities from objects in class: ",class(object),call.=FALSE)
  else
  p
}

predictSurvOrder.mfp <- function(object,newdata,times,...){
  # require(mfp)
  p <- predictSurvOrder.coxph(object$fit,newdata=newdata,times=times)
  p
}


predictSurvOrder.survnnet <- function(object,newdata,times,train.data,...){
  require(Design)
  learndat <- train.data
  learndat$nnetFactor <- predict(object,train.data,...)
  newdata$nnetFactor <- predict(object,newdata)
  nnet.form <- reformulate("nnetFactor",object$call$formula[[2]])
  fit.nnet <- cph(nnet.form,data=learndat,se.fit=FALSE,surv=TRUE,x=TRUE,y=TRUE)
  p <- predictSurvOrder.cph(fit.nnet,newdata=newdata,times=times)
  p
}


#predictSurvOrder.rpart <- function(object,newdata,times,train.data,...){
#  # require(rpart)
#  require(Design)
#  learndat <- train.data
#  nclass <- length(unique(object$where))
#  learndat$rpartFactor <- factor(predict(object,newdata=train.data,...))
#  newdata$rpartFactor <- factor(predict(object,newdata=newdata))
#  rpart.form <- reformulate("rpartFactor",object$call$formula[[2]])
#  #  fit.rpart <- cph(rpart.form,data=learndat,se.fit=FALSE,surv=TRUE,x=TRUE,y=TRUE)
#  fit.rpart <- prodlim(rpart.form,data=learndat)
#  p <- predictSurvOrder(fit.rpart,newdata=newdata,times=times)
#  p
#}


predictSurvOrder.coxph <- function(object,newdata,times,...){
  require(survival)
  p <- c(predict(object,newdata=newdata))
  #  p <- predictSurvProb(object,newdata=newdata,times=times)[,1]
  if (NROW(p) != NROW(newdata))
    stop("Prediction failed")
  # order(-p)
  -p
}


predictSurvOrder.cph <- function(object,newdata,times,...){
  p <- predict(object,newdata=newdata)
  if (NROW(p) != NROW(newdata))
    stop("Prediction failed")
  #    order(-p)
  -p
}

#predictSurvOrder.prodlim <- function(object,newdata,times,...){
#  require(prodlim)
#  p <- predict(object=object,type="surv",newdata=newdata,times=times,mode="matrix",level.chaos=1)
#  p[is.na(p)] <- 0
#  if (is.null(dim(p)))
#    {if (length(p)!=length(times))
#       stop("Prediction failed")}
#  else{
#    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
#      stop("Prediction failed")
#  }
#  p
#}


## library randomSurvivalForest
## predictSurvOrder.rsf <- function(object,newdata,times,...){
##   p <- predict.rsf(object,newdata=newdata,times=times,bytimes=TRUE,fill="last")
##   if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
##     stop("Prediction failed")
##   p
## }


predictSurvOrder.psm <- function(object,newdata,times,...){
  p <- predict(object,newdata=newdata)
  if (NROW(p) != NROW(newdata))
    stop("Prediction failed")
  p
}

#predictSurvOrder.phnnet <- function(object,newdata,times,train.data,...){
##  require(survnnet)
#  learndat <- train.data
#  seeds <- sample(1:1000,size=10)
#  object$call$data <- learndat
#  re.fitter <- lapply(seeds,function(s){
#    set.seed(s)
#    refit <- eval(object$call)
#    list(learn=predict(refit,learndat),
#         val=predict(refit,newdata))
#  })
#  learndat$nnetFactor <- rowMeans(do.call("cbind",lapply(re.fitter,function(x)x[["learn"]])))
#  newdata$nnetFactor <- rowMeans(do.call("cbind",lapply(re.fitter,function(x)x[["val"]])))
#  #  a <- predict(object,train.data,...)
#  #  b <- predict(object,newdata)
#  #  print(cbind(do.call("cbind",lapply(re.fitter,function(x)x[["learn"]]))[1:10,],learndat$nnetFactor[1:10]))
#  #  stop()
#  #  learndat$nnetFactor <- predict(object,train.data,...)
#  #  newdata$nnetFactor <- predict(object,newdata)
#  nnet.form <- reformulate("nnetFactor",object$call$formula[[2]])
#  require(Design)
#  fit.nnet <- cph(nnet.form,data=learndat,se.fit=FALSE,surv=TRUE,x=TRUE,y=TRUE)
#  p <- predictSurvOrder.cph(fit.nnet,newdata=newdata,times=times)
#  p
#}


