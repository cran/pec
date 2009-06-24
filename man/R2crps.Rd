\name{R2}
\alias{R2}
\alias{crps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing prediction error curves}
\description{
 Cumulating prediction error curves over time and
 a time-dependent $R^2$ like measure.  
}
\usage{
R2(object, who, what, times ,nullModel=1)
crps(object, who, what, times, start)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object with estimated prediction error curves
    obtained with the function \link{pec}}
  \item{who}{Which models in \code{object$models} should be
    considered. Default is all models for crps and all models except the
  reference for R2}
  \item{what}{Usually \code{pred.error} -- if crossvalidation or bootstrap methods are used
    then also the apparent.error and if \code{replan=boot.632plus} then
    also the BootB0.error, the estimated overfit, and the  NoInf.error
    can be cumulated or compared to a reference prediction model.}
  \item{times}{Time points at which the summaries are shown.}
  \item{start}{Only for crps: the time point at which the cumulation is started}
  \item{nullModel}{
    Position of the model whose prediction error is used
    as the reference in the denominator when
    constructing $R^2$}
}
\details{
  The cumulative prediction error (continuous ranked probability score)
  is defined as the area under the prediction error curve.

  In survival analysis the prediction error of the Kaplan-Meier
  estimator plays a similar role as the total sum of squares in linear regression.
  Hence, it is a sensible reference model for $R^2$.  }

\value{A matrix with a column for every requested prediction model}
% \references{ ~put references to the literature/web site
%   here ~ }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}

\seealso{\code{\link{pec}}}

\examples{set.seed(18713)

\dontrun{
dat=SimSurv(100)
nullmodel=prodlim(Hist(time,status)~1,data=dat)
pmodel=coxph(Surv(time,status)~X1+X2,data=dat)
perror=pec(list(KaplanMeier=nullmodel,Cox=pmodel),Hist(time,status)~1,data=dat)

## cumulative prediction error
crps(perror,times=1) # between min time and 1
crps(perror,times=1,start=0) # between 0 and 1
crps(perror,times=seq(0,1,.2),start=0) # between 0 and seq(0,1,.2)

R2(perror,times=seq(0,1,.1))
}
}

% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

\keyword{survival}% at least one, from doc/KEYWORDS
