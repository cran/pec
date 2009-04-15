roc.Surv <- function(x,y,times,bandwidth="smooth"){
  time <- y[,"time"]
  status <- y[,"status"]
  survorder <- order(time, - status)
  time <- time[survorder]
  status <- status[survorder]
  x <- x[survorder]
  N <- length(x)
  ##   test 
  ##   bandwidth <- N^(-1/3)
  ##   test
  if (missing(times)) times <- sort(unique(time))
  breax <- sort(unique(x))
  # conditional survival probability
  ## condSurv.x <- prodlim(Surv(time,status)~NN(x),bandwidth=bandwidth)
  condSurv.x <- prodlim(Surv(time,status)~x,bandwidth=bandwidth)
  ##   summary(condSurv.x)
  Surv.TX <- predict(condSurv.x,newdata=data.frame(x=x),times=times,type="surv",mode="matrix",bytime=TRUE,level.chaos=0)
  sortx <- sort(x)
  tmpX <- do.call("rbind",lapply(breax,function(c){sortx > c})) # matrix with elements 0 and 1
  rownames(tmpX) <- paste("c=",breax,sep="")
  colnames(tmpX) <- paste("Xi=",sort(x),sep="")
  ## biSurv <- lapply(1:length(times),function(s){apply(Surv.TX[s,]*tmpX,1,sum)/N})
  Surv.T <- mean(condSurv.x,newdata= data.frame(x=x),times=times)$surv
  Surv.marg <- predict(prodlim(Surv(time,status)~1),type="surv",times=times)
  Surv.X <- predict(prodlim(Surv(x,rep(1,N))~1),times=breax,level.chaos=0,type="surv")  ## empirical marker dist
  names(Surv.X) <- paste("X=",breax,sep="")
##   print(Surv.T)
  ##   print(Surv.marg)
  out <- lapply(1:length(times),function(s){
    ##     print(list(s=s))
    surv.tx.s <- Surv.TX[s,]
    surv.tx.s[is.na(surv.tx.s)] <- 0
    ##     print(t(1*(tmpX)))
    biSurv <- rowMeans(t(surv.tx.s*t(tmpX)),na.rm=TRUE)
    ## Sens <- (Surv.X - biSurv)/(1-Surv.marg[s])
    Sens <- (Surv.X - biSurv)/(1-mean(surv.tx.s))
##     if (s==2){
##       print(Surv.X)
##       print(biSurv)
##       print(Surv.T[s])
##     }
    Spec <- 1-biSurv/Surv.T[s]
    Auc <- auc.default(Sens,Spec)
    list("Sensitivity"=Sens,"Specificity"=Spec,"Auc"=Auc)
  })
  names(out) <- paste("t",times,sep="=")
  out
}

## y <- 1:10
## s <- rep(1,10)
## x <- c(2,4,3,1,5,8,10,9,6,7)
## r <- roc.Surv(x,Surv(y,s),times=5)

## R <- roc(x,y>5)


## r[3]
## y.t3 <- y>3
## cbind(roc.default(x,y.t3)$Sens,r[[3]]$Sens)

## x <- roc.Surv(Surv(dfs.time,rep(1,157))~weight,data=epo)
## names(x)[131]
## stat.131 <- epo$dfs.time>54.7
## roc.131 <- roc.default(epo$weight,stat.131)
## cbind(roc.131$Sens,x[[131]]$Sens)
