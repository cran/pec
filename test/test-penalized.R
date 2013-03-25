library(pec)
library(penalized)
source("predictSurvProb.penfitS3.R")
data(nki70)

# -------------------------Fitting the model-------------------------

## S4 fit
pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
                 unpenalized = ~ER+Age+Diam+N+Grade, data = nki70, lambda1 = 1)

## my way of calling yields the same thing packed in an S3 object:
penS3 <- penalizedS3(Surv(time,event)~ER+Age+Diam+pen(8:77)+N+Grade,
                     data=nki70,
                     lambda1=1)
## or
penS3 <- penalizedS3(Surv(time,event)~ER+Age+Diam+pen(TSPYL5,Contig63649_RC)+pen(10:77)+N+Grade,
                     data=nki70,
                     lambda1=1)
## also this works
penS3 <- penalizedS3(Surv(time,event)~ER+Age+pen(8:33)+Diam+pen(34:77)+N+Grade,
                     data=nki70,
                     lambda1=1)


## pre-optimized

library(penalized)
library(pec)
data(nki70)
penLopt <- penalizedOpt(Surv(time, event)~ER+Age+Diam+N+Grade + pen(8:77),data = nki70,optL1.trace=FALSE,optL1.minlambda1=0,optL1.minlambda1=15,optL1.model="cox",optL1.fold=50,optL1.standardize=TRUE,penalized.standardize=TRUE)
penL1 <- penalizedS3(Surv(time, event)~pen(8:77) + ER+Age+Diam+N+Grade, data = nki70,lambda1=1)
cox <- coxph(Surv(time, event)~ER+Age+Diam+N+Grade, data = nki70)
apparrentC <- cindex(list(optLambda1=penLopt,fixedLambda1=penL1,cox=cox),formula=Surv(time, event)~ER+Age+Diam+N+Grade,data=nki70)
## need to set B up and adjust M to your data set
set.seed(17)
internalValC <- cindex(list(optLambda1=penLopt,fixedLambda1=penL1,cox=cox),formula=Surv(time, event)~ER+Age+Diam+N+Grade,data=nki70,splitMethod="bootcv",B=3,M=round(.632*NROW(nki70)))

# ------------------------predicting survival------------------------

predictSurvProb.penfitS3(penS3,newdata=nki70[43:48,],times=c(0,13,18,22))
# ----------------------prediction error curves----------------------

apparentPec <- pec(list(penalized=penS3),
                formula=Surv(time,event)~1,
                data=nki70)
plot(apparent)

# bootstrap-crossvalidation (subsampling training size M=100)
set.seed(13)
bootcvPec <- pec(list(penalized=penS3),
              formula=Surv(time,event)~1,
              data=nki70,
              splitMethod="bootcv",
              B=10,
              M=100)
plot(bootcv)

# ----------------------cindex ----------------------------------

## a single evaluation time

apparentC13 <- cindex(list(penalized=penS3),
                    formula=Surv(time,event)~1,
                    data=nki70,eval.times=13)

## discrimination curve over time
apparentC <- cindex(list(penalized=penS3),
                    formula=Surv(time,event)~1,
                    data=nki70,eval.times=1:13)
plot(apparentC)

# bootstrap-crossvalidation (subsampling training size M=100)
set.seed(13)
bootcvC <- cindex(list(penalized=penS3),
                    formula=Surv(time,event)~1,
                    data=nki70,splitMethod="bootcv",B=10,M=100,eval.times=1:13)
plot(bootcvC)
