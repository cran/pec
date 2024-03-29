# {{{ roxy header
#' Estimation of censoring probabilities
#' 
#' This function is used internally by the function \code{pec} to obtain
#' inverse of the probability of censoring weights.
#' 
#' Inverse of the probability of censoring weights (IPCW) usually refer to the
#' probabilities of not being censored at certain time points. These
#' probabilities are also the values of the conditional survival function of
#' the censoring time given covariates. The function ipcw estimates the
#' conditional survival function of the censoring times and derives the
#' weights.
#' 
#' IMPORTANT: the data set should be ordered, \code{order(time,-status)} in
#' order to get the values \code{IPCW.subjectTimes} in the right order for some
#' choices of \code{method}.
#' 
#' @aliases ipcw ipcw.none ipcw.marginal ipcw.nonpar ipcw.cox ipcw.aalen
#' @param formula A survival formula like, \code{Surv(time,status)~1}, where
#' as usual status=0 means censored. The status variable is internally
#' reversed for estimation of censoring rather than survival
#' probabilities. Some of the available models (see argument
#' \code{model}) will use predictors on the right hand side of the
#' formula.
#' @param data The data used for fitting the censoring model
#' @param method Censoring model used for estimation of the
#' (conditional) censoring distribution.
#' @param args A list of arguments which is passed to method
#' @param times For \code{what="IPCW.times"} a vector of times at
#' which to compute the probabilities of not being censored.
#' @param subjectTimes For \code{what="IPCW.subjectTimes"} a vector of
#' individual times at which the probabilities of not being censored
#' are computed.
#' @param subjectTimesLag If equal to \code{1} then obtain
#' \code{G(T_i-|X_i)}, if equal to \code{0} estimate the conditional
#' censoring distribution at the subjectTimes,
#' i.e. (\code{G(T_i|X_i)}).
#' @param what Decide about what to do: If equal to
#' \code{"IPCW.times"} then weights are estimated at given
#' \code{times}.  If equal to \code{"IPCW.subjectTimes"} then weights
#' are estimated at individual \code{subjectTimes}.  If missing then
#' produce both.
#' @return \item{times}{The times at which weights are estimated}
#' \item{IPCW.times}{Estimated weights at \code{times}}
#' \item{IPCW.subjectTimes}{Estimated weights at individual time values
#' \code{subjectTimes}} \item{fit}{The fitted censoring model}
#' \item{method}{The method for modelling the censoring distribution}
#' \item{call}{The call}
#' @author Thomas A. Gerds \email{tag@@biostat.ku.dk}
#' @seealso \code{\link{pec}}
#' @keywords survival
#' @examples
#' 
#' library(prodlim)
#' library(rms)
#' library(survival)
#' dat=SimSurv(30)
#' 
#' dat <- dat[order(dat$time),]
#' 
#' # using the marginal Kaplan-Meier for the censoring times
#' 
#' WKM=ipcw(Hist(time,status)~X2,
#'   data=dat,
#'   method="marginal",
#'   times=sort(unique(dat$time)),
#'   subjectTimes=dat$time)
#' plot(WKM$fit)
#' WKM$fit
#' 
#' # using the Cox model for the censoring times given X2
#' library(survival)
#' WCox=ipcw(Hist(time=time,event=status)~X2,
#'   data=dat,
#'   method="cox",
#'   times=sort(unique(dat$time)),
#'   subjectTimes=dat$time)
#' WCox$fit
#' 
#' plot(WKM$fit)
#' lines(sort(unique(dat$time)),
#'       1-WCox$IPCW.times[1,],
#'       type="l",
#'       col=2,
#'       lty=3,
#'       lwd=3)
#' lines(sort(unique(dat$time)),
#'       1-WCox$IPCW.times[5,],
#'       type="l",
#'       col=3,
#'       lty=3,
#'       lwd=3)
#' 
#' # using the stratified Kaplan-Meier
#' # for the censoring times given X2
#' 
#' WKM2=ipcw(Hist(time,status)~X2,
#'   data=dat,
#'   method="nonpar",
#'   times=sort(unique(dat$time)),
#'   subjectTimes=dat$time)
#' plot(WKM2$fit,add=FALSE)
#' 
#'
# }}}
# {{{ method ipcw
#' @export
ipcw <- function(formula,
                 data,
                 method,
                 args,
                 times,
                 subjectTimes,
                 subjectTimesLag=1,
                 what){
    if (!missing(what))
        stopifnot(all(match(what,c("IPCW.times","IPCW.subjectTimes"))))
    if (missing(what) || match("IPCW.times",what,nomatch=FALSE)){
        stopifnot(length(times)>0)
    }
    class(method) <- method
    UseMethod("ipcw",method)
}
# }}}
# {{{ uncensored data: return just 1

##' @export
ipcw.none <- function(formula,data,method,args,times,subjectTimes,subjectTimesLag,what){
    if (missing(subjectTimesLag)) subjectTimesLag=1
    if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
    call <- match.call()
    #  weigths at requested times
    if (match("IPCW.times",what,nomatch=FALSE)){
        IPCW.times <- rep(1,length(times))
    }
    else
        IPCW.times <- NULL
    #  weigths at subject specific event times
    if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
        IPCW.subjectTimes <- rep(1,length(subjectTimes))
    }
    else
        IPCW.subjectTimes <- NULL
    out <- list(times=times,IPCW.times=IPCW.times,IPCW.subjectTimes=IPCW.subjectTimes,call=call,method=method)
    class(out) <- "IPCW"
    out
}

# }}}
# {{{ reverse Random Survival Forests
##' @export
ipcw.rfsrc <- function(formula,data,method,args,times,subjectTimes,subjectTimesLag,what){
    if (!(requireNamespace("randomForestSRC",quietly=TRUE)))
        stop("Need package randomForestSRC for this, but is not available here.")
    if (missing(subjectTimesLag)) subjectTimesLag=1
    if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
    call <- match.call() ## needed for refit in crossvalidation loop
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       specials=NULL,
                                       unspecialsDesign=FALSE)
    wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
    ## wdata <- as.data.frame(EHF)
    wdata$status <- 1-wdata$status
    wform <- update(formula,"Surv(time,status)~.")
    ## require(randomForestSRC)
    stopifnot(NROW(na.omit(wdata))>0)
    if (missing(args) || is.null(args))
        ## args <- list(bootstrap="none",ntree=1000,nodesize=NROW(data)/2)
        args <- list(ntree=1000)
    ## if (is.null(args$importance) & (args$importance!="none"))
    args$importance <- "none"
    fit <- do.call(randomForestSRC::rfsrc,c(list(wform,data=wdata),args))
    ## print(fit)
    fit$call <- NULL
    #  weigths at requested times
    #  predicted survival probabilities for all training subjects are in object$survival
    #  out-of-bag prediction in object$survival.oob
    if (match("IPCW.times",what,nomatch=FALSE)){
        IPCW.times <- predictSurvProb(fit,newdata=wdata,times=times)
    }
    else
        IPCW.times <- NULL
    #  weigths at subject specific event times
    if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
        pmat <- fit$survival
        jtimes <- fit$time.interest
        IPCW.subjectTimes <- sapply(1:length(subjectTimes),function(i){
            Ci <- subjectTimes[i]
            pos <- prodlim::sindex(jump.times=jtimes,eval.times=Ci,comp="smaller",strict=(subjectTimesLag==1))
            c(1,pmat[i,])[1+pos]
        })
    }
    else
        IPCW.subjectTimes <- NULL
    out <- list(times=times,
                IPCW.times=IPCW.times,
                IPCW.subjectTimes=IPCW.subjectTimes,
                fit=fit,
                call=call,
                method=method)
    ## browser()
    ## print(head(IPCW.subjectTimes))
    class(out) <- "IPCW"
    out
}
##' @export
ipcw.forest <- function(formula,data,method,args,times,subjectTimes,subjectTimesLag,what){
    if (!(requireNamespace("randomForestSRC",quietly=TRUE)))
        stop("Need package randomForestSRC for this, but is not available here.")
    if (missing(subjectTimesLag)) subjectTimesLag=1
    if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
    call <- match.call() ## needed for refit in crossvalidation loop
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       specials=NULL,
                                       unspecialsDesign=FALSE)
    wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
    ## wdata$status <- 1-wdata$status
    wform <- update(formula,"Surv(time,status)~.")
    ## require(randomForestSRC)
    stopifnot(NROW(na.omit(wdata))>0)
    if (missing(args) || is.null(args))
        args <- list(ntree=1000)
    args$importance <- "none"
    fit <- do.call(randomForestSRC::rfsrc,c(list(wform,data=wdata),args))
    ## print(fit)
    fit$call <- NULL
    # forest weights
    FW <- predict(fit,newdata=wdata,forest.wt=TRUE)$forest.wt
    #  weigths at requested times
    #  predicted survival probabilities for all training subjects are in object$survival
    #  out-of-bag prediction in object$survival.oob
    if (match("IPCW.times",what,nomatch=FALSE)){
        # reverse Kaplan-Meier with forest weigths
        IPCW.times <- apply(data,1,function(i){
            predict(prodlim::prodlim(Hist(time,status)~1,data=wdata,reverse=TRUE,caseweights=FW[i,]),times=times)
        })
    }
    else
        IPCW.times <- NULL
    #  weigths at subject specific event times
    if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
        IPCW.subjectTimes <- sapply(1:length(subjectTimes),function(i){
            ## browser()
            prodlim::predictSurvIndividual(prodlim::prodlim(Hist(time,status)~1,data=wdata,reverse=TRUE,caseweights=FW[i,]),lag=1)[i]
        })
    }
    else
        IPCW.subjectTimes <- NULL
    out <- list(times=times,
                IPCW.times=IPCW.times,
                IPCW.subjectTimes=IPCW.subjectTimes,
                fit=fit,
                call=call,
                method=method)
    class(out) <- "IPCW"
    out
}
# }}}
# {{{ reverse Kaplan-Meier
##' @export
ipcw.marginal <- function(formula,data,method,args,times,subjectTimes,subjectTimesLag,what){
    if (missing(subjectTimesLag)) subjectTimesLag=1
    if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
    call <- match.call()
    formula <- update.formula(formula,"~1")
    fit <- prodlim::prodlim(formula,data=data,reverse=TRUE)
    #  weigths at requested times
    if (match("IPCW.times",what,nomatch=FALSE)){
        IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
    }
    else
        IPCW.times <- NULL
    #  weigths at subject specific event times
    if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
        IPCW.subjectTimes <- prodlim::predictSurvIndividual(fit,lag=subjectTimesLag)
    }
    else
        IPCW.subjectTimes <- NULL
    out <- list(times=times,
                IPCW.times=IPCW.times,
                IPCW.subjectTimes=IPCW.subjectTimes,
                fit=fit,
                call=call,
                method=method)
    class(out) <- "IPCW"
    out
    ##   locsubjectTimes <- match(subjectTimes,fit$time,nomatch=NA)
    ##   if (any(is.na(locsubjectTimes))) stop("Can not locate all individual observation times" )
    ##   IPCW.subjectTimes <- c(1,fit$surv)[locsubjectTimes] ## at (subjectTimes_i-)
    ##   IPCW.times <- c(1,fit$surv)[prodlim::sindex(jump.times=fit$time,eval.times=times) +1] ## at all requested times
}
# }}}
# {{{ reverse Stone-Beran 
##' @export
ipcw.nonpar <- function(formula,data,method,args,times,subjectTimes,subjectTimesLag,what){
    if (missing(subjectTimesLag)) subjectTimesLag=1
    if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
    call <- match.call()
    fit <- prodlim::prodlim(formula,data=data,reverse=TRUE,bandwidth="smooth")
    #  weigths at requested times
    if (match("IPCW.times",what,nomatch=FALSE)){
        IPCW.times <- predict(fit,newdata=data,times=times,level.chaos=1,mode="matrix",type="surv")
    }
    else
        IPCW.times <- NULL
    #  weigths at subject specific event times
    if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
        IPCW.subjectTimes <- prodlim::predictSurvIndividual(fit,lag=subjectTimesLag)
    }
    else
        IPCW.subjectTimes <- NULL
    out <- list(times=times,
                IPCW.times=IPCW.times,
                IPCW.subjectTimes=IPCW.subjectTimes,
                fit=fit,
                call=call,
                method=method)
    class(out) <- "IPCW"
    out
}
# }}}
# {{{ reverse Cox via Harrel's package
##' @export
ipcw.cox <- function(formula,data,method,args,times,subjectTimes,subjectTimesLag,what){
    ## require(rms)
    if (missing(subjectTimesLag)) subjectTimesLag=1
    if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
    call <- match.call()
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       specials=c("strat"),
                                       stripSpecials=c("strat"),
                                       specialsDesign=FALSE,
                                       unspecialsDesign=FALSE)
    if (is.null(EHF$strat))
        wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
    else
        wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design,EHF$strat))
    ## wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
    wdata$status <- 1-wdata$status
    wform <- update(formula,"Surv(time,status)~.")
    stopifnot(NROW(na.omit(wdata))>0)
    if (missing(args) || is.null(args))
        args <- list(x=TRUE,y=TRUE,eps=0.000001)
    args$surv <- TRUE
    fit <- do.call(rms::cph,c(list(wform,data=wdata),args))
    ## fit <- rms::cph(wform,data=wdata,surv=TRUE,x=TRUE,y=TRUE)
    #  weigths at requested times
    if (match("IPCW.times",what,nomatch=FALSE)){
        IPCW.times <- rms::survest(fit,newdata=wdata,times=times,se.fit=FALSE)$surv
    }
    else
        IPCW.times <- NULL
    #  weigths at subject specific event times
    if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
        if (subjectTimesLag==1)
            IPCW.subjectTimes <- rms::survest(fit,times=subjectTimes-min(diff(c(0,unique(subjectTimes))))/2,what='parallel')
        else if (subjectTimesLag==0){
            IPCW.subjectTimes <- rms::survest(fit,times=subjectTimes,what='parallel')
        }
        else stop("SubjectTimesLag must be 0 or 1")
    }
    else
        IPCW.subjectTimes <- NULL
    out <- list(times=times,
                IPCW.times=IPCW.times,
                IPCW.subjectTimes=IPCW.subjectTimes,
                fit=fit,
                call=call,
                method=method)
    class(out) <- "IPCW"
    out
}
# }}}
# {{{ reverse Aalen method via the timereg package
##' @export
ipcw.aalen <- function(formula,data,method,args,times,subjectTimes,subjectTimesLag,what){
    if (missing(subjectTimesLag)) subjectTimesLag=1
    if (missing(what)) what=c("IPCW.times","IPCW.subjectTimes")
    call <- match.call()
    EHF <- prodlim::EventHistory.frame(formula,
                                       data,
                                       specials=NULL,
                                       unspecialsDesign=FALSE)
    wdata <- data.frame(cbind(unclass(EHF$event.history),EHF$design))
    ## wdata <- as.data.frame(EHF)
    wdata$status <- 1-wdata$status
    wform <- update(formula,"Surv(time,status)~.")
    stopifnot(NROW(na.omit(wdata))>0)
    fit <- do.call(timereg::aalen,list(formula=formula,data=wdata,n.sim=0))
    fit$call <- NULL
    #  weigths at requested times
    if (match("IPCW.times",what,nomatch=FALSE)){
        IPCW.times <- predictSurvProb(fit,newdata=wdata,times=times)
    }  else {
        IPCW.times <- NULL
    }
    if (match("IPCW.subjectTimes",what,nomatch=FALSE)){
        if (subjectTimesLag==1) 
            IPCW.subjectTimes <- diag(predictSurvProb(fit,newdata=data,times=pmax(0,subjectTimes-min(diff(unique(subjectTimes)))/2)))
        else if (subjectTimesLag==0)
            IPCW.subjectTimes <- diag(predictSurvProb(fit,newdata=data,times=subjectTimes))
        else stop("SubjectTimesLag must be 0 or 1")
    }
    else
        IPCW.subjectTimes <- NULL
    out <- list(times=times,IPCW.times=IPCW.times,IPCW.subjectTimes=IPCW.subjectTimes,fit=fit,call=call,method=method)
    class(out) <- "IPCW"
    out
}
# }}}

