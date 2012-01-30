ccr <- function(time,status,risk){
  N <- length(time)
  pairs <- data.frame(do.call("rbind",lapply(1:N,function(i){
    do.call("rbind",lapply(1:N,function(j){
      c(time[i],time[j],status[i],status[j],risk[i],risk[j],risk[i]>risk[j],time[i]<time[j])
    }))
  })))
  names(pairs) <- c("time.i","time.j","D.i","D.j","r.i","r.j","riskOrdered","timeOrdered")
  pairs
}

library(pec)
library(riskRegression)
library(cmprsk)
d <- data.frame(time=1:4,status=c(1,2,0,1))
risk <- matrix(c(2,1,3,4))
## pairs:
## id1 vs id2 in, concordant: Y[1]<Y[2] & D[1]=1
## id1 vs id3 in, discordant: risk[1]<risk[3]
## id2 vs id3 out: D[2]=2
## id2 vs id1 out: D[2]=2
## id3 vs id1 out: Y[3]>Y[1] & D[1]=1
## id3 vs id2 in, concordant: Y[3]>=Y[2] & D[2]=2

cindex(risk,formula=Hist(time,status)~1,data=d)
concordance.CR(risk,d$time,d$status,failcode=1,cencode=0)
concordance.simple.CR(risk,d$time,d$status,failcode=1,cencode=0)

x <-  cindex(list(risk,risk1,risk2),formula=Hist(time,status)~1,data=d)

N <- 4
dd <- prodlim:::SimCompRisk(N)
r1 <- matrix(sample(1:N))
r2 <- matrix(sample(N:1))

ccr(dd$time,dd$cause,r1)

dd <- data.frame(time=c(3,2,4,1),status=c(1,1,1,1),cause=c(2,1,1,1),r1=c(1,2,3,4))
r1 <- matrix(c(1,3,2,4))
sdd <- dd[order(dd$time),]

cindex(matrix(dd$r1),formula=Hist(time,cause)~1,data=dd)
cindex(matrix(sdd$r1),formula=Hist(time,cause)~1,data=sdd)

concordance.CR(sdd$r1,sdd$time,sdd$cause,failcode=1,cencode=0)
concordance.CR(dd$r1,dd$time,dd$cause,failcode=1,cencode=0)

concordance.simple.CR(dd$r1,dd$time,dd$cause,failcode=1,cencode=0)

data(Melanoma)
## Melanoma <- Melanoma[40:80,]
## Melanoma <- Melanoma[Melanoma$status!=0,]
f <- FGR(Hist(time,status)~age+thick+ulcer,data=Melanoma)
f2 <- FGR(Hist(time,status)~age+ulcer,data=Melanoma)
f3 <- FGR(Hist(time,status)~age,data=Melanoma)
cindex(list(f,f2,f3),formula=Hist(time,status)~1,data=Melanoma)

cindex(list(f,f2,f3),formula=Hist(time,status)~1,data=Melanoma,splitMethod="bootcv",B=10)

cindex(list(f,f2,f3),formula=Hist(time,status)~age+thick+ulcer,data=Melanoma,cens.model="cox",splitMethod="bootcv",B=10)

library(riskRegression)
library(cmprsk)
d <- get(load(file="~/tmp/train.b"))
f <- FGR(Hist(time,status)~age+thick+ulcer,data=d)

p <- predictEventProb(f,newdata=Melanoma,times=1500)
concordance.CR(p,Melanoma$time,Melanoma$status,failcode=1,cencode=0)
concordance.simple.CR(p,Melanoma$time,Melanoma$status,failcode=1,cencode=0)
u <- cindex(list(f),formula=Hist(time,status)~1,data=Melanoma)

u <- cindex(list(f),formula=Hist(time,status)~1,data=Melanoma,splitMethod="bootcv",B=3,M=NROW(Melanoma))
u <- cindex(list(f),formula=Hist(time,status)~1,data=Melanoma,splitMethod="bootcv",B=3,M=NROW(Melanoma),verbose=FALSE)

M1 <- Melanoma[Melanoma$status!=0,]
p <- predictEventProb(f,newdata=M1,times=1500)
cindex(list(f),formula=Hist(time,status)~1,data=M1)
concordance.CR(p,M1$time,M1$status,failcode=1,cencode=0)
concordance.simple.CR(p,M1$time,M1$status,failcode=1,cencode=0)


cindex(list(f2,f),formula=Hist(time,status)~1,data=Melanoma)
cindex(list(f,f2),formula=Hist(time,status)~1,data=Melanoma)
cindex(list(f,f2,f3),formula=Hist(time,status)~1,data=Melanoma,eval.time=1500)
