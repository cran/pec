# {{{ RSF GBSG2 reproducible example
library("pec")
library("rms")
library("randomSurvivalForest")
library("party")
data(GBSG2)
GBSG2$status <- GBSG2$cens
GBSG2 <- GBSG2[order(GBSG2$time,-GBSG2$status),]
GBSG2$grade.bin <- as.factor(as.numeric(GBSG2$tgrade!="I"))
levels(GBSG2$grade.bin) <- c("I","II/III")
fitformGBSG=Surv(time,status)~age+tsize+pnodes+progrec+estrec+grade.bin
fitcox=cph(fitformGBSG,data=GBSG2,surv=TRUE,se.fit=FALSE)
set.seed(17)
fitrsf=rsf(fitformGBSG,data=GBSG2,forest=TRUE,ntree=100)
set.seed(17)
fitcforest <- pecCforest(fitformGBSG, data=GBSG2,
controls=cforest_classical(ntree=100))
extends <- function(...)TRUE

set.seed(2006)
fitpec <- pec(list("Cox"=fitcox,"rsf"=fitrsf,"cforest"=fitcforest),formula=Surv(time,cens)~age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2,cens.model="cox",splitMethod="Boot632plus",maxtime=2000,B=5,keep.index=TRUE,keep.matrix=TRUE, verbose=TRUE)
## set.seed(2006)
## fitpec2 <- pec(list("Cox"=fitcox,"rsf"=fitrsf,"cforest"=fitcforest),formula=Surv(time,cens)~age+tsize+grade.bin+pnodes+progrec+estrec,data=GBSG2,cens.model="cox",splitMethod="Boot632plus",maxtime=2000,B=5,keep.index=TRUE,keep.matrix=TRUE, verbose=TRUE,doMC=TRUE)
crps.t2000 <- crps(fitpec,times=2000)
crps.t2000

            ## AppErr BootCvErr NoInfErr Boot632plusErr
## KaplanMeier  0.180     0.182    0.180          0.180
## Cox          0.157     0.160    0.197          0.159
## rsf          0.106     0.157    0.209          0.146
## cforest      0.143     0.159    0.202          0.155

            ## AppErr BootCvErr NoInfErr Boot632plusErr
## KaplanMeier  0.180     0.182    0.180          0.180
## Cox          0.157     0.160    0.197          0.159
## rsf          0.106     0.158    0.209          0.146
## cforest      0.143     0.160    0.202          0.155

# }}}
