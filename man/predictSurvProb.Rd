\name{predictSurvProb}
\alias{predictSurvProb}
\alias{predictSurvProb.aalen}
\alias{predictSurvProb.cox.aalen}
\alias{predictSurvProb.coxph}
\alias{predictSurvProb.cph}
\alias{predictSurvProb.default}
\alias{predictSurvProb.glm}
\alias{predictSurvProb.matrix}
\alias{predictSurvProb.mfp}
\alias{predictSurvProb.prodlim}
\alias{predictSurvProb.psm}
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
\method{predictSurvProb}{cox.aalen}(object, newdata, times,...)
\method{predictSurvProb}{coxph}(object, newdata, times,...)
\method{predictSurvProb}{cph}(object, newdata, times,...)
\method{predictSurvProb}{default}(object, newdata, times,...)
\method{predictSurvProb}{glm}(object, newdata, times,...)
\method{predictSurvProb}{matrix}(object, newdata, times,...)
\method{predictSurvProb}{mfp}(object, newdata, times,...)
\method{predictSurvProb}{prodlim}(object, newdata, times,...)
\method{predictSurvProb}{psm}(object, newdata, times,...)
\method{predictSurvProb}{survfit}(object, newdata, times,...)
\method{predictSurvProb}{phnnet}(object, newdata, times,train.data,...)
\method{predictSurvProb}{survnnet}(object, newdata, times,train.data,...)
\method{predictSurvProb}{rpart}(object, newdata, times,train.data,...)

}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{A model for which predicted probabilities are desired.}
  \item{newdata}{A data frame containing predictor variable combinations for
    which predictions are desired}
  \item{times}{A vector of times in the range of the response
    variable, e.g. times when the response is a survival object, at
    which the exceedance probabilities (i.e. the survival probabilities) are
    returned.}
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
  Design), coxph (library survival), aalen (library timereg),
  cox.aalen (library timereg), mfp (library mfp), phnnet (library
  survnnet), survnnet (library survnnet), rpart (library rpart),
  product.limit (library prodlim), survfit (library survival), psm
  (library Design), glm (library stats).
}

\value{
  A matrix with as many rows as \code{NROW(newdata)} and as many columns
  as \code{length(times)}. Each entry should be a probability and in
  rows the values should be decreasing.}
% \references{ ~put references to the literature/web site here ~ }

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
% \examples{

% }

\keyword{survival}% at least one, from doc/KEYWORDS
