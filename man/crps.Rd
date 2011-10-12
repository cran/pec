\name{crps}
\alias{crps}
\alias{ibs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Summarizing prediction error curves}
\description{
  Computes the cumulative prediction error curves, aka integrated Brier
  scores, in ranges of time.
}
\usage{
crps(object, models, what, times, start)
ibs(object, models, what, times, start)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{An object with estimated prediction error curves
    obtained with the function \link{pec}}
  \item{models}{Which models in \code{object$models} should be
    considered.
  }
  \item{what}{
    The name of the entry in \code{x} to be cumulated. Defauls to \code{PredErr} 
    Other choices are \code{AppErr}, \code{BootCvErr}, \code{Boot632},
    \code{Boot632plus}.
  }
  \item{times}{Time points at which the integration of the prediction
    error curve stops.}
  \item{start}{The time point at which the integration of the prediction
    error curve is started.}
}
\details{
  The cumulative prediction error (continuous ranked probability score)
  is defined as the area under the prediction error curve, hence the
  alias name, ibs, which is short for integrated Brier score. 

}
\value{A matrix with a column for the crps (ibs) at every requested time point and a row
  for each model}
\references{
    E. Graf et al.  (1999),
  Assessment and comparison of prognostic classification schemes for
  survival data. Statistics in Medicine, vol 18, pp= 2529--2545.

  Gerds TA, Cai T & Schumacher M (2008)
  The performance of risk prediction models
  Biometrical Journal, 50(4), 457--479

}
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}

\seealso{\code{\link{pec}}}

\examples{set.seed(18713)

dat=prodlim:::SimSurv(100)
nullmodel=prodlim(Hist(time,status)~1,data=dat)
pmodel=coxph(Surv(time,status)~X1+X2,data=dat)
perror=pec(list(KaplanMeier=nullmodel,Cox=pmodel),Hist(time,status)~1,data=dat)

## cumulative prediction error
crps(perror,times=1) # between min time and 1
## same thing:
ibs(perror,times=1) # between min time and 1
crps(perror,times=1,start=0) # between 0 and 1
crps(perror,times=seq(0,1,.2),start=0) # between 0 and seq(0,1,.2)

}

\keyword{survival}% at least one, from doc/KEYWORDS
