library(pec)
set.seed(130971)

# ------------------------simulate some data------------------------

dat <- SimSurv(100)

# fit some candidate Cox models and compute the Kaplan-Meier estimate

Models <- list("Cox.X1"=coxph(Surv(time,status)~X1,data=dat),
               "Cox.X2"=coxph(Surv(time,status)~X2,data=dat),
               "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=dat))

f <- coxph(Surv(time,status)~X2,data=dat)
predictSurvProb(f,newdata=data.frame(X2=0:1),times=c(0,150,200,245))

k <- prodlim(Surv(time,status)~1,data=dat)
predictSurvProb(k,newdata=data.frame(X2=0:1),times=c(0,150,200,245,300))


# compute the prediction error with the old method
set.seed(127)
system.time(PredError <- pec(object=Models,formula=Surv(time,status)~X1+X2,data=dat,exact=TRUE,cens.model="marginal",replan="bootcv",B=10,M=66,verbose=TRUE))

# compute the prediction error with the new method and do the multiSplitTest
check.code("pec")
set.seed(127)
out <- pecStar.list(object=Models[1],
                    formula=Surv(time,status)~X1+X2,
                    data=dat,
                    exact=TRUE,
                    cens.model="marginal",
                    replan="bootcv",
                    B=10,
                    M=66,
                    verbose=TRUE,
                    keepResiduals=TRUE,
                    keep.index=TRUE,
                    multiSplitTest=TRUE,
                    testTimes=c(10,50,200),
                    testIBS=c(0,60))

## ttt <- out$splitMethod$index[,1]
## ddd <- dat[ttt,]
## vvv <- dat[!((1:NROW(dat)) %in% ttt),]
## kmkmkm <- prodlim(Surv(time,status)~1,data=ddd)
## predictSurvProb(kmkmkm,times=c(10,20,30,200,3000),newdata=vvv)
## fff <- coxph(Surv(time,status)~X1,data=ddd)
## round(predictSurvProb(fff,times=c(10,20,30,200),newdata=vvv),2)


system.time(PredErrorStar <- pecStar.list(object=Models,
                                          formula=Surv(time,status)~X1+X2,
                                          data=dat,
                                          exact=TRUE,
                                          cens.model="marginal",
                                          replan="bootcv",
                                          B=10,
                                          M=66,
                                          verbose=TRUE,
                                          multiSplitTest=TRUE,
                                          multiSplitTest.times=30))

##                   ____     __        ___      _ _____         _   
## __   ____ _ _ __ |  _ \  __\ \      / (_) ___| |_   _|__  ___| |_ 
## \ \ / / _` | '_ \| | | |/ _ \ \ /\ / /| |/ _ \ | | |/ _ \/ __| __|
##  \ V / (_| | | | | |_| |  __/\ V  V / | |  __/ | | |  __/\__ \ |_ 
##   \_/ \__,_|_| |_|____/ \___| \_/\_/  |_|\___|_| |_|\___||___/\__|

# median p-values of 10 bootstrap cross-validation steps

PredErrorStar$multiSplitTest$multiSplitTest.times
plot(PredErrorStar,xlim=c(0,100))
t(t(PredErrorStar$multiSplitTest$pVal))
hist(PredErrorStar$multiSplitTest$pairedPvalues[1,])
hist(PredErrorStar$multiSplitTest$pairedPvalues[2,])
hist(PredErrorStar$multiSplitTest$pairedPvalues[3,])

# the comparison Cox X1 versus Cox X2 is not significant because
# the alternative is "Cox X2 has a smaller prediction error than Cox X1


# -------------------------------checks-------------------------------

# compute the prediction error with the new method and do NOT do multiSplitTest
set.seed(127)
system.time(PredErrorStar2 <- pec:::pecStar.list(object=Models,formula=Surv(time,status)~X1+X2,data=dat,exact=TRUE,cens.model="marginal",replan="bootcv",B=10,M=66,verbose=TRUE,multiSplitTest=FALSE))


# check if all apparent errors are the same
sapply(1:4,function(m){
  p2=PredErrorStar$AppErr[[m]]
  p1=PredError$AppErr[[m]]
  all(p1[!is.na(p1)]==p2[!is.na(p2)])
})

# check if all apparent errors are the same
sapply(1:4,function(m){
  p2=PredErrorStar2$AppErr[[m]]
  p1=PredError$AppErr[[m]]
  all(p1[!is.na(p1)]==p2[!is.na(p2)])
})

# check if all the bootstrapCrossValidation  errors are the same
sapply(1:4,function(m){
  p2=PredErrorStar$PredErr[[m]]
  p1=PredError$PredErr[[m]]
  all(round(p1[!is.na(p1)],12)==round(p2[!is.na(p2)],12))
})

# check if all the bootstrapCrossValidation  errors are the same
sapply(1:4,function(m){
  p2=PredErrorStar2$PredErr[[m]]
  p1=PredError$PredErr[[m]]
  all(round(p1[!is.na(p1)],12)==round(p2[!is.na(p2)],12))
})
