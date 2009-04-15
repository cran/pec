\name{cindex}
\alias{cindex}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Computes the concordance index for right censored survival time data}
\description{
  The concordance index (C-index) is the frequency of concordant pairs
  among all usuable pairs, where a pair consists of the response and a
  continuous marker of one individual. Without censoring the statistic
  is related to the Wilcoxon statistic and Somers coefficient. The
  C-index equals the area under the ROC curve when the response is binary.
  This function provides the inverse of the probability of censoring
  weigthing estimate to adjust for right censoring . 
}
\usage{
cindex(object,
       data,model
       cens.formula = formula,
       cens.model = c("cox", "marginal", "nonpar", "aalen", "none"),
       conf.int = NULL,
       B = 1000,
       na.action=na.omit,
       tau)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A named list of prediction models. Each entry
    is either an R-object for which a \link{predictSurvOrder}
    method exists (see details) or
    a \code{call} that evaluates to such an R-object, see examples.
    For resampling (e.g. when \code{replan} is "boot632plus")
    all the models in this list must include their
    \code{call} in their value. 
  }
  \item{data}{a data frame that contains all the variables in
    \code{formula} and \code{cens.formula}}
  \item{cens.formula}{a formula of the form \code{Surv(time,status) ~
      z1+z2+...} where z1 and z2 are confounders.}
  \item{cens.model}{Method for estimating inverse probability of
    censoring weigths: \code{cens.model="nonpar"} means nonparametric
    kernel smoothing of the Kaplan-Meier estimator for the censoring
    variable using symmetric nearest neighborhoods -- this is possible
    for arbitrary many strata variables but at most one continuous
    covariate, see the function prodlim in the package prodlim. A regression
    when (\code{cens.model="aalen"}) or (\code{cens.model="cox"}).
    \code{marginal}: ignoring all covariates Kaplan-Meier estimator 
  }
  \item{conf.int}{level of confidence: a decimal between 0 and 1, i.e.\ .95}
  \item{B}{Number of bootstrap replications}

  \item{na.action}{
    Passed to \code{\link{model.frame}}
  }
  
  \item{tau}{The upper limit of integration for estimates of restricted
    C-index. see details section}
}
\details{
  The C-index can be represented by four integrals: two for the event
  time and two for the continuous marker. The range of the event time
  variable may be restricted by specifying the upper limit of
  integration tau. 
 
}
\value{
The C-index 
}
\references{Harrel et al. 1996, Statistics in Medicine, vol 15, pages 361-87}
\author{Thomas A Gerds \email{tag@biostat.ku.dk}}

% \examples{
% library(survival)
% dat=SimSurv(300)
% cox=coxph(Surv(time,status)~X1+X2,data=dat)
% cox$call
% cindex(list("Cox"=cox),formula=Surv(time,status)~X1+X2,data=dat)
% }

\keyword{survival}
