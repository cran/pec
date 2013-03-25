# {{{ testing seed
library(pec)
set.seed(17)
d <- SimSurv(100)
f <- coxph(Surv(time,status)~X2,data=d)
set.seed(13)
a=pec(f,splitMethod="bootcv",B=3,M=63,keep.index=TRUE,verbose=F)
b=a$splitMethod$index
set.seed(13)
c=pec:::resolvesplitMethod(splitMethod="bootcv",N=100,M=63,B=3)$index
stopifnot(all.equal(b,c))
# }}}

## library(party)
## h <- cforest(Surv(time,status)~X2,data=d)
## isS4(h)
## f <- cph(Surv(time,status)~X2,data=d,surv=TRUE)
## set.seed(19)
## A <- pec(f,splitMethod="cv5",B=1,M=63,keep.index=TRUE,verbose=F)

# {{{ testing ipcw
set.seed(18)
A=pec(f,B=30,splitMethod="bootcv")
set.seed(18)
A1=pec(f,B=30,ipcw.refit=T,splitMethod="bootcv")
cbind(A$BootCvErr$CoxModel,A1$BootCvErr$CoxModel)
## plot(A,xlim=c(0,100))
## plot(A1,add=TRUE,lty=3)
B=pec(f,cens.model="cox")

# }}}

# {{{ testing splitMethods: cvk, loocv

# Comparing the predictive performance of some standard
# survival regression models 
# --------------------------------------------------------------------

library(pec)
library(rms) # thanks to Frank Harrell 
data(GBSG2)
GBSG2$status <- GBSG2$cens
GBSG2 <- GBSG2[order(GBSG2$time,-GBSG2$status),]
GBSG2$grade.bin <- as.factor(as.numeric(GBSG2$tgrade!="I"))
levels(GBSG2$grade.bin) <- c("I","II/III")

GBSG2=GBSG2[sample(1:NROW(GBSG2),size=50),]
# Kaplan-Meier (survival package)
system.time(pbc.fit0 <- survfit(Surv(time,status)~1,data=GBSG2))
# Kaplan-Meier (prodlim package) is faster
system.time(pbc.fit0 <- prodlim(Hist(time,status)~1,data=pbc) )
## Cox model (rms)
Cox=cph(Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2,x=TRUE,y=TRUE,surv=TRUE,se.fit=FALSE)

set.seed(17)
check.code("pec")
loocv <- pec.list(object=list(Cox),formula=Surv(time,status)~1,data=GBSG2,cens.model="marginal",splitMethod="loocv",B=1,keep.matrix=TRUE,verbose=TRUE)
cv5 <- pec.list(object=list(Cox),formula=Surv(time,status)~1,data=GBSG2,cens.model="marginal",splitMethod="cv5",B=2,keep.matrix=TRUE,verbose=TRUE)
# }}}

# {{{ testing splitMethods: boot632, boot632+

## d <- SimSurv(300)
## f2 <- coxph(Surv(time,status)~rcs(X1)+X2,data=d)

library(pec)
library(rms) # thanks to Frank Harrell 
data(GBSG2)
GBSG2$status <- GBSG2$cens
GBSG2 <- GBSG2[order(GBSG2$time,-GBSG2$status),]
GBSG2$grade.bin <- as.factor(as.numeric(GBSG2$tgrade!="I"))
levels(GBSG2$grade.bin) <- c("I","II/III")

d <- SimSurv(100)
library(randomForestSRC)

RSF=rfsrc(Surv(time,status)~age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2,forest=TRUE)
save(RSF,file="~/research/SoftWare/pec/exa/RSF-gbsg2.rda")
set.seed(17)
b632 <- pec.list(object=list(RSF),formula=Surv(time,status)~1,data=GBSG2,cens.model="marginal",splitMethod="boot632plus",B=10,M=500)
set.seed(17)
b632a <- pec.list(object=list(RSF),noinf.permute=T,formula=Surv(time,status)~1,data=GBSG2,cens.model="marginal",splitMethod="boot632plus",B=10,M=500)
plot(b632a$NoInfErr[[2]],b632$NoInfErr[[2]])

f3 <- rfsrc(Survrsf(time,status)~X1+X2,data=d)
set.seed(17)
b632 <- pec.list(object=list(f3),formula=Surv(time,status)~1,data=d,cens.model="marginal",splitMethod="boot632plus",B=10)
set.seed(17)
b632a <- pec.list(object=list(f3),noinf.permute=T,formula=Surv(time,status)~1,data=d,cens.model="marginal",splitMethod="boot632plus",B=10)
plot(b632a$NoInfErr[[1]],b632$NoInfErr[[1]])
plot(b632a$NoInfErr[[2]],b632$NoInfErr[[2]])
## plot(b632,xlim=c(0,100))
## plot(b632a,add=T,lty=3)
# }}}



         
