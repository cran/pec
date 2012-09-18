library(pec)
data(cost)
f <- coxph(Surv(time,status)~age+sex,data=cost)
fitpec <- pec(list("cox" = f),data = cost,formula = Surv(time, status) ~ 1,splitMethod = "Boot632plus",B = 11,M = 350,keep.index = TRUE,keep.matrix = TRUE)
## fitpec
## Range of integration: 0 and time=4259 :
## Integrated Brier score (crps):
##                IBS[0;time=4259]
## KaplanMeier            0.200
## cox                    0.177
stopifnot(round(summary(fitpec$AppErr[[1]]),4)==c(0,0.1716,0.1889,0.1861,0.2352,0.25))
stopifnot(round(summary(fitpec$AppErr[[2]]),4)==c(0,0.1362,0.17,0.1634,0.208,0.2223))
stopifnot(round(summary(fitpec$Boot632plusErr[[1]]),4)==c(0,0.1706,0.1874,0.1849,0.2346,0.25,1))
stopifnot(round(summary(fitpec$Boot632plusErr[[2]]),4)==c(0,0.1394,0.1697,0.1641,0.2078,0.2246,1))


