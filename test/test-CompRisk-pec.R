library(pec)
library(riskRegression)
library(cmprsk)
data(Melanoma)
Melanoma$logthick <- log(Melanoma$thick)
arr.fit <- ARR(Hist(time,status)~sex+epicel+ulcer+age+logthick,data=Melanoma,cause=1,cens.model="cox",cens.formula=~sex+epicel+ulcer+age+logthick)
cox.fit <- CSC(list(cause1=Hist(time,status)~sex+epicel+ulcer+age+logthick,cause2=Hist(time,status)~sex+age),data=Melanoma)
fg.fit <- FGR(Hist(time,status)~sex+epicel+ulcer+age+logthick,data=Melanoma,cause=1)
lrr.fit <- LRR(Hist(time,status)~sex+epicel+ulcer+age+logthick,data=Melanoma,cause=1,cens.model="cox",cens.formula=~sex+epicel+ulcer+age+logthick)

mlist <- list(AbsRisk=arr.fit,CauseSpecCox=cox.fit,FineGray=fg.fit,LogisticRisk=lrr.fit)
perr.link <- pec(mlist,formula=Hist(time,status)~sex+epicel+ulcer+age+logthick,cens.model="cox",data=Melanoma,maxtime=2500)
plot(perr.link)

perr.link1 <- pec(mlist,formula=Hist(time,status)~1,times=1499,exact=FALSE,start=NULL,cause=1,cens.model="marginal",data=Melanoma[1:10,],maxtime=2500,reference=FALSE)
perr.link2 <- pec(mlist,formula=Hist(time,status)~1,times=c(1499,1500),exact=FALSE,start=NULL,cause=1,cens.model="marginal",data=Melanoma[1:10,],maxtime=2500,reference=FALSE)
cbind(perr.link1$AppErr,sapply(perr.link2$AppErr,function(x)x[1]))


library(survival)
library(riskRegression)
library(cmprsk)
library(pec)
data(pbc)
f1  <- CSC(Hist(time,status)~sex+edema,data=pbc)
## f1  <- CSC(Hist(time,status)~sex+edema,cause=2,data=pbc)
f2  <- CSC(Hist(time,status)~sex,data=pbc)
f3  <- FGR(Hist(time,status)~sex+edema,cause=2,data=pbc)
f4  <- FGR(Hist(time,status)~sex+edema,cause=2,data=pbc)
p1 <- pec(list(f1,f2,f3,f4),formula=Hist(time,status)~1,data=pbc,cause=2)

predictEventProb(f1,newdata=pbc[17,],times=1000)
predictEventProb(f1,newdata=pbc[17,],times=1000,cause=2)
predictEventProb(f1,newdata=pbc[17,],times=1000,cause=1)



perr.link1 <- pec(list(LogisticRisk=lrr.fit),formula=Hist(time,status)~sex+epicel+ulcer+age+logthick,times=1499,exact=FALSE,start=NULL,cause=1,cens.model="cox",data=Melanoma,maxtime=2500)
perr.link2 <- pec(list(LogisticRisk=lrr.fit),formula=Hist(time,status)~sex+epicel+ulcer+age+logthick,times=c(1499,1500),exact=FALSE,start=NULL,cause=1,cens.model="cox",data=Melanoma,maxtime=2500)
