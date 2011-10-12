library(pec)
library(penalized)
source("predictSurvProb.penfitS3.R")
data(nki70)

# -------------------------Fitting the model-------------------------

## S4 fit
pen <- penalized(Surv(time, event), penalized = nki70[,8:77],
                 unpenalized = ~ER+Age+Diam+N+Grade, data = nki70, lambda1 = 10)

## my way of calling yields the same thing packed in an S3 object:
penS3 <- penalizedS3(Surv(time,event)~ER+Age+Diam+pen(8:77)+N+Grade,
                     data=nki70,
                     lambda1=1)
## or
penS3 <- penalizedS3(Surv(time,event)~ER+Age+Diam+pen(TSPYL5,Contig63649_RC)+pen(10:77)+N+Grade,
                     data=nki70,
                     lambda1=1)
## even this works
penS3 <- penalizedS3(Surv(time,event)~ER+Age+pen(8:33)+Diam+pen(34:77)+N+Grade,
                     data=nki70,
                     lambda1=1)

# ------------------------predicting survival------------------------

predictSurvProb.penfitS3(penS3,newdata=nki70[43:48,],times=c(0,13,18,22))

apparent <- pec(list(penalized=penS3),
                formula=Surv(time,event)~1,
                data=nki70)
plot(apparent)
set.seed(13)
bootcv <- pec(list(penalized=penS3),
              formula=Surv(time,event)~1,
              data=nki70,
              replan="bootcv",
              B=10,
              M=100)
plot(bootcv)


