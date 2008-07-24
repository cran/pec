\name{ipcw}
\alias{ipcw}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimation of censoring probabilities}
\description{This function is used internally by the function \code{pec}
  to obtain inverse of the probability of censoring weights.}
\usage{
ipcw(formula,
     data,
     model = c("cox", "marginal", "nonpar", "aalen", "none"),
     times,
     otimes)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{A survival formula like, Surv(time,status)~1,
    where as usual status=0 means
    censored. The status variable is internally reversed for estimation
    of censoring rather than survival probabilities. Some of the
    available models (see argument \code{model})
  will use predictors on the right hand side of the formula.}
  \item{data}{The data used for fitting the censoring model}
  \item{model}{Censoring model used for estimation of the (conditional) censoring
    distribution.}
  \item{times}{Time values at which the probabilities of not being
    censored are returned}
  \item{otimes}{The time values, usually the original event times in the data.frame \code{data},
    at which lagged individual probabilities of not being censored are
    returned. 
  }
}
\details{Inverse of the probability of censoring weights (IPCW) usually 
  refer to the probabilities of not being censored at certain time
  points. These probabilities are also the values of the conditional survival
  function of the censoring time given covariates. The function ipcw
  estimates the conditional survival
  function of the censoring times and derives the weights.
  
  Currently the data set should be order by time to get the value \code{wt.obs} in
  the right order for some choices of \code{model}!
}
\value{
  \item{wt}{Estimated weights at \code{times}}
  \item{wt.obs}{Estimated weights at individual time values \code{otimes}}
  \item{fit}{The fitted censoring model}
}
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}

\seealso{\code{\link{pec}}}
\examples{
library(survival)
data(pbc)

pbc <- pbc[order(pbc$time),]

# using the marginal Kaplan-Meier for the censoring times

WKM=ipcw(Surv(time,status)~edema,data=pbc,model="marginal",times=sort(unique(pbc$time)),otimes=pbc$time)
plot(WKM$fit)
WKM$fit

# using the Cox model for the censoring times given edema

WCox=ipcw(Surv(time,status)~edema,data=pbc,model="cox",times=sort(unique(pbc$time)),otimes=pbc$time)
WCox$fit

plot(WKM$fit)
lines(sort(unique(pbc$time)),1-WCox$wt[1,],type="l",col=2,lty=3,lwd=3)
lines(sort(unique(pbc$time)),1-WCox$wt[2,],type="l",col=3,lty=3,lwd=3)

# using the stratified Kaplan-Meier for the censoring times given edema

WKM2=ipcw(Surv(time,status)~edema,data=pbc,model="nonpar",times=sort(unique(pbc$time)),otimes=pbc$time)
plot(WKM2$fit,add=TRUE)

 }
\keyword{survival}% at least one, from doc/KEYWORDS
