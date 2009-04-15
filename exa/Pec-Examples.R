# Comparing the predictive performance of some standard
# survival regression models 
# --------------------------------------------------------------------

library(pec)
library(Design) # thanks Frank Harrell 

## ,----
## |  Predicting the time to death for patients with
## |  primary biliary cirrhosis
## `----

data(pbc)

## ,----
## | 0. Benchmark prediction:
## | Kaplan-Meier predicts the same probabilities for all patients
## | 1. The bilirubin level is a predictor?
## | 2. The log(bilirubin level) predicts better?
## | 3. Michael WK always uses splines
## | 4. A multifactor model predicts always better!
## |    covariate selection according to book Therneau & Grambsch
## | 5. Michael WK always uses splines 
## `----

# Kaplan-Meier (survival package)
pbc.fit0 <- survfit(Surv(time,status)~1,data=pbc)

# Kaplan-Meier (prodlim package) is faster
pbc.fit0 <- prodlim(Hist(time,status)~1,data=pbc) 

tmp <- coxph(Surv(time,status)~age+edema+bili+protime+alb,data=pbc)

see <- pec.list(object=list(fullcox=tmp),
                formula=Surv(time,status)~age+edema+bili+protime+alb,
                data=pbc,
                maxtime=3000,
                exact=F,
                cens.model="marginal",
                replan="boot632plus",
                B=100,
                keep.matrix=TRUE,
                na.accept=10,
                verbose=TRUE)



# Cox regression models
pbc.fit1 <- cph(Surv(time,status)~bili,data=pbc,surv=T,x=T,y=T) 
pbc.fit2 <- cph(Surv(time,status)~log(bili),data=pbc,surv=T)
pbc.fit3 <- cph(Surv(time,status)~rcs(bili,parms=3),data=pbc,surv=T)
pbc.fit4 <- cph(Surv(time,status)~age+edema+log(bili)+log(protime)+log(alb),data=pbc,surv=T)
pbc.fit5 <- cph(Surv(time,status)~age+edema+rcs(bili,parms=3)+rcs(protime,parms=3)+rcs(alb,parms=3),data=pbc,surv=T)

## Evaluate the predictive performance with the Brier score

pbc.models <- list("Kaplan-Meier"=pbc.fit0,
                   "Cox bili"=pbc.fit1,
                   "Cox log(bili)"=pbc.fit2,
                   "Cox rcs(bili,3)"=pbc.fit3,
                   "Cox multi"=pbc.fit4,
                   "Cox rcs(multi,3)"=pbc.fit5)
## ,----
## |  Dependent censoring!!!
## |  To avoid bias we find censoring weights 
## |  conditional on the covariates.
## `----

set.seed(20092007)
pbc.pec <- pec.list(object=pbc.models,
                    formula=Surv(time,status)~age+edema+bili+protime+alb,
                    data=pbc,
                    maxtime=3000,
                    exact=TRUE,
                    cens.model="cox",
                    replan="boot632plus",
                    B=100,
                    keep.matrix=TRUE,
                    na.accept=10,
                    verbose=TRUE)

print(pbc.pec)
summary(pbc.pec,times=c(300,600,900,1200,1500,1800))

plot(pbc.pec)
plot(pbc.pec,who=c(1,2,4,6),smooth=TRUE)
plot(pbc.pec,who=c(1,2,4,6),smooth=TRUE,what="AppErr")
plot(pbc.pec,who=c(1,2,4,6),smooth=TRUE,what="OutOfBagErr")
plot(pbc.pec,who=c(1,2,4,6),smooth=TRUE,what="special",specials=list(what=c("PredErr"),bench=1))


library(timereg)

         
