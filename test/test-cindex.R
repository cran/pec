# {{{ load libraries
library(pec)
library(Hmisc)
library(survival)
# }}}
# {{{ generate dummy data
set.seed(16)
dat <- SimSurv(30,cens=FALSE)
datC <- SimSurv(30,cens=TRUE)
# }}}
# {{{ fit dummy models
fit12 <- coxph(Surv(time,status)~X1+X2,data=datC)
fit1 <- coxph(Surv(time,status)~X1,data=datC)
fit2 <- coxph(Surv(time,status)~X2,data=datC)
# }}}
# {{{ compare C to Harrell's C uncensored data
Cpec <- cindex(list("Cox X1+X2"=fit12,"Cox X1"=fit1,"Cox X2"=fit2),formula=Surv(time,status)~1,data=dat,eval.times=Inf)
p1 <- predictSurvProb(fit1,newdata=dat,times=100)
p2 <- predictSurvProb(fit2,newdata=dat,times=100)
p12 <- predictSurvProb(fit12,newdata=dat,times=100)
harrellC1 <- rcorr.cens(p1,with(dat,Surv(time,status)))
harrellC2 <- rcorr.cens(p2,with(dat,Surv(time,status)))
harrellC12 <- rcorr.cens(p12,with(dat,Surv(time,status)))
stopifnot(harrellC1[["C Index"]]==Cpec$AppCindex[["Cox X1"]])
stopifnot(harrellC1[["Relevant Pairs"]]==Cpec$Pairs[["Cox X1"]])
stopifnot(harrellC1[["Concordant"]]==Cpec$Concordant[["Cox X1"]])
stopifnot(harrellC2[["C Index"]]==Cpec$AppCindex[["Cox X2"]])
stopifnot(harrellC2[["Relevant Pairs"]]==Cpec$Pairs[["Cox X2"]])
stopifnot(harrellC2[["Concordant"]]==Cpec$Concordant[["Cox X2"]])
stopifnot(harrellC12[["C Index"]]==Cpec$AppCindex[["Cox X1+X2"]])
stopifnot(harrellC12[["Relevant Pairs"]]==Cpec$Pairs[["Cox X1+X2"]])
stopifnot(harrellC12[["Concordant"]]==Cpec$Concordant[["Cox X1+X2"]])
message("DONE:cindex for uncensored data equal to Hmisc rcorr.cens")
##
## tied predictions in and out
##
dummy <- data.frame(time=1:4,status=rep(1,4),x=c(1,2,1,3))
## u <- unlist(lapply(1:4,function(i)lapply(1:4,function(j)if (i!=j) cbind(dummy[i,],dummy[j,]) else NULL)),rec=F)
## u <- do.call("rbind",u[!sapply(u,is.null)])
harrellC3 <- rcorr.cens(dummy$x,with(dummy,Surv(time,status)),outx=FALSE)
harrellC3a <- rcorr.cens(dummy$x,with(dummy,Surv(time,status)),outx=TRUE)
C3 <- cindex(list("x"=dummy$x),formula=Surv(time,status)~1,data=dummy,eval.times=Inf,tiedPredictionsIn=TRUE)
C3a <- cindex(list("x"=dummy$x),formula=Surv(time,status)~1,data=dummy,eval.times=Inf,tiedPredictionsIn=FALSE)
stopifnot(C3$AppCindex[["x"]]==harrellC3["C Index"])
stopifnot(C3a$AppCindex[["x"]]==harrellC3a["C Index"])
message("DONE:test effect of ties in predictions on c-index ")
##
## tied outcome in and out
##
dummy2 <- data.frame(time=c(1,1,2,2,3,3,4,4,5,5),status=rep(1,10),x=c(1,2,1,3,2,4,3,5,4,5),z=c(1,2.1,1.2,3,2,4.1,3.1,5.1,4,5))
harrellC4a <- rcorr.cens(dummy2$x,with(dummy2,Surv(time,status)),outx=TRUE)
C4a <- cindex(list("x"=dummy2$x),formula=Surv(time,status)~1,data=dummy2,eval.times=Inf,tiedPredictionsIn=FALSE,tiedOutcomeIn=TRUE)
C4b <- cindex(list("x"=dummy2$x),formula=Surv(time,status)~1,data=dummy2,eval.times=Inf,tiedPredictionsIn=FALSE,tiedOutcomeIn=FALSE)
stopifnot(C4a$AppCindex[["x"]]!=harrellC4a["C Index"])
stopifnot(C4b$AppCindex[["x"]]==harrellC4a["C Index"])
message("DONE:test effect of ties in outcome on c-index ")
# }}}
# {{{ compare C to Harrell's C with censored data
#
# Harrells index ignores pairs where Y[i]==Y[j]
# this seems inappropriate if 
# a) status[i]=1 & status[j]=0 (since then uncensored T[j] > T[i])
# or
# b) pred[i] == pred[j] and status[i]==status[j]==1 (since then pred and outcome are concordant)
#
dummy3 <- data.frame(time=c(1,1,2),status=c(1,1,1),x=c(2,2.2,3))
harrellC5 <- rcorr.cens(dummy3$x,with(dummy3,Surv(time,status)),outx=FALSE)
C5a <- cindex(list("x"=dummy3$x),formula=Surv(time,status)~1,data=dummy3,eval.times=Inf,tiedPredictionsIn=TRUE,tiedOutcomeIn=TRUE)
C5b <- cindex(list("x"=dummy3$x),formula=Surv(time,status)~1,data=dummy3,eval.times=Inf,tiedPredictionsIn=TRUE,tiedOutcomeIn=FALSE)
# }}}
