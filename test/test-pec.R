library(pec)
data(cost)
f <- coxph(Surv(time,status)~age+sex,data=cost)
set.seed(17)
fitpec <- pec(list("cox" = f),data = cost,formula = Surv(time, status) ~ 1,splitMethod = "Boot632plus",B = 11,M = 350,keep.index = TRUE,keep.matrix = TRUE)
## fitpec
## Range of integration: 0 and time=4259 :
## Integrated Brier score (crps):
##                IBS[0;time=4259]
## KaplanMeier            0.200
## cox                    0.177
stopifnot(round(summary(fitpec$AppErr[[1]]),4)==c(0,0.1716,0.1889,0.1861,0.2352,0.2500))
stopifnot(round(summary(fitpec$AppErr[[2]]),4)==c(0,0.1362,0.17,0.1634,0.208,0.2223))
stopifnot(round(summary(fitpec$Boot632plusErr[[1]]),4)==c(0,0.1659,0.1891,0.1855,0.2346,0.25,2))
stopifnot(round(summary(fitpec$Boot632plusErr[[2]]),4)==c(0,0.1321,0.1712,0.1641,0.2096,0.2254,2))

## f <- coxph(Surv(time,status)~age+sex,data=cost)
## x <- pec(list("Cox"=f),formula=Surv(time,status)~1,data=cost,B=1,splitMethod="cv10",keep.index=TRUE,keep.matrix=TRUE)
## x <- pec(list("Cox"=f),formula=Surv(time,status)~1,data=cost,B=10,splitMethod="bootcv",keep.index=TRUE,keep.matrix=TRUE)
## f1 <- cph(Surv(time,status)~age+sex,data=cost,surv=TRUE)
## x <- pec(list("Cox"=f1),formula=Surv(time,status)~1,data=cost,B=10,splitMethod="cv10",keep.index=TRUE,keep.matrix=TRUE)


