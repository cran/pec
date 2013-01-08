\name{predictSurvProb}
\alias{predictSurvProb}
\alias{predictSurvProb.aalen}
\alias{predictSurvProb.riskRegression}
\alias{predictSurvProb.cox.aalen}
\alias{predictSurvProb.coxph}
\alias{predictSurvProb.cph}
\alias{predictSurvProb.default}
\alias{predictSurvProb.rsf}
\alias{predictSurvProb.matrix}
\alias{predictSurvProb.mfp}
\alias{predictSurvProb.pecCforest}
\alias{predictSurvProb.prodlim}
\alias{predictSurvProb.psm}
\alias{predictSurvProb.selectCox}
\alias{predictSurvProb.survfit}
\alias{predictSurvProb.survnnet}
\alias{predictSurvProb.phnnet}
\alias{predictSurvProb.rpart}


%- Also NEED an '\alias' for EACH other topic documented here.
\title{Predicting survival probabilities}
\description{
  Function to extract survival probability predictions
  from various modeling
  approaches. The most prominent one is the Cox regression model which can be
  fitted for example with `coxph' and with `cph'.
  
  The function predictSurvProb is a generic function that means it invokes
  specifically designed functions depending on the 'class' of the first argument.
}
\usage{
 predictSurvProb(object, newdata, times, ...)
\method{predictSurvProb}{aalen}(object, newdata, times,...)
\method{predictSurvProb}{riskRegression}(object, newdata, times,...)
\method{predictSurvProb}{cox.aalen}(object, newdata, times,...)
\method{predictSurvProb}{coxph}(object, newdata, times,...)
\method{predictSurvProb}{cph}(object, newdata, times,...)
\method{predictSurvProb}{default}(object, newdata, times,...)
\method{predictSurvProb}{rsf}(object, newdata, times,...)
\method{predictSurvProb}{matrix}(object, newdata, times,...)
\method{predictSurvProb}{mfp}(object, newdata, times,...)
\method{predictSurvProb}{pecCforest}(object, newdata, times,...)
\method{predictSurvProb}{prodlim}(object, newdata, times,...)
\method{predictSurvProb}{psm}(object, newdata, times,...)
\method{predictSurvProb}{selectCox}(object, newdata, times,...)
\method{predictSurvProb}{survfit}(object, newdata, times,...)
\method{predictSurvProb}{phnnet}(object, newdata, times,train.data,...)
\method{predictSurvProb}{survnnet}(object, newdata, times,train.data,...)
\method{predictSurvProb}{rpart}(object, newdata, times,train.data,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{A fitted model from which to extract predicted survival probabilities}
 \item{newdata}{A data frame containing predictor variable combinations for
   which to compute predicted survival probabilities.
 }
  \item{times}{A vector of times in the range of the response
    variable, e.g. times when the response is a survival object, at
    which to return the survival probabilities.
  }
   \item{train.data}{An optional data frame which contains the
     response and predictor variable combinations
     in which the prediction model was trained}
  \item{\dots}{Additional arguments that are passed on to the current method.
  }
}
\details{
  The function \code{pec} requires survival probabilities for each row in newdata
  at requested times. These probabilities are extracted from
  a fitted model of class \code{CLASS} with the function
  \code{predictSurvProb.CLASS}.
  
  Currently there are \code{predictSurvProb} methods for objects of class cph (library
  rms), coxph (library survival), aalen (library timereg),
  cox.aalen (library timereg), mfp (library mfp), phnnet (library
  survnnet), survnnet (library survnnet), rpart (library rpart),
  product.limit (library prodlim), survfit (library survival), psm
  (library rms)
}

\value{
  A matrix with as many rows as \code{NROW(newdata)} and as many columns
  as \code{length(times)}. Each entry should be a probability and in
  rows the values should be decreasing.}
\references{
  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012). Evaluating Random Forests
  for Survival Analysis Using Prediction Error Curves. Journal of Statistical Software,
  50(11), 1-23. URL http://www.jstatsoft.org/v50/i11/.
}

\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}
\note{
  In order to assess the predictive performance of a new survival model a
  specific \code{predictSurvProb} S3 method has to be written. For
  examples, see the bodies of the
  existing methods.

  The performance of the assessment procedure,
  in particular for resampling 
  where the model is repeatedly evaluated, will be improved by
  supressing in the call to the model all the computations that are not
  needed for probability prediction. For example, \code{se.fit=FALSE}
  can be set in the call to \code{cph}.  
}

\seealso{\code{\link{predict}},\code{\link{survfit}}}
\examples{
# generate some survival data
library(prodlim)
d <- SimSurv(100)
# then fit a Cox model
library(rms)
coxmodel <- cph(Surv(time,status)~X1+X2,data=d,surv=TRUE)

# predicted survival probabilities can be extracted
# at selected time-points:
ttt <- quantile(d$time)
# for selected predictor values:
ndat <- data.frame(X1=c(0.25,0.25,-0.05,0.05),X2=c(0,1,0,1))
# as follows
predictSurvProb(coxmodel,newdata=ndat,times=ttt)

## simulate some learning and some validation data
learndat <- SimSurv(100)
valdat <- SimSurv(100)
## use the learning data to fit a Cox model
fitCox <- coxph(Surv(time,status)~X1+X2,data=learndat)
## suppose we want to predict the survival probabilities for all patients
## in the validation data at the following time points:
## 0, 12, 24, 36, 48, 60
psurv <- predictSurvProb(fitCox,newdata=valdat,times=seq(0,60,12))
## This is a matrix with survival probabilities
## one column for each of the 5 time points
## one row for each validation set individual

\dontrun{
# the same can be done e.g. for a randomSurvivalForest model
library(randomSurvivalForest)
rsfmodel <- rsf(Survrsf(time,status)~X1+X2,data=d)
predictSurvProb(rsfmodel,newdata=ndat,times=ttt)
}
}

\keyword{survival}% at least one, from doc/KEYWORDS
