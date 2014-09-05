\name{R2}
\alias{R2}
\title{Explained variation for survival models}
\usage{
R2(object, models, what, times, reference = 1)
}
\arguments{
  \item{object}{An object with estimated prediction error
  curves obtained with the function \link{pec}}

  \item{models}{For which of the models in
  \code{object$models} should we compute $R^2(t). By
  default all models are used except for the reference
  model.}

  \item{what}{The name of the entry in \code{x} to be used.
  Defauls to \code{PredErr} Other choices are
  \code{AppErr}, \code{BootCvErr}, \code{Boot632},
  \code{Boot632plus}.}

  \item{times}{Time points at which the summaries are
  shown.}

  \item{reference}{Position of the model whose prediction
  error is used as the reference in the denominator when
  constructing $R^2$}
}
\value{
A matrix where the first column holds the times and the
following columns are the corresponding $R^2$ values for
the requested prediction models.
}
\description{
This function computes a time-dependent $R^2$ like measure
of the variation explained by a survival prediction model,
by dividing the mean squared error (Brier score) of the
model by the mean squared error (Brier score) of a
reference model which ignores all the covariates.
}
\details{
In survival analysis the prediction error of the
Kaplan-Meier estimator plays a similar role as the total
sum of squares in linear regression. Hence, it is a
sensible reference model for $R^2$.
}
\examples{
set.seed(18713)
library(prodlim)
library(survival)
dat=SimSurv(100)
nullmodel=prodlim(Hist(time,status)~1,data=dat)
pmodel1=coxph(Surv(time,status)~X1+X2,data=dat)
pmodel2=coxph(Surv(time,status)~X2,data=dat)
perror=pec(list(Cox1=pmodel1,Cox2=pmodel2),Hist(time,status)~1,data=dat,reference=TRUE)
R2(perror,times=seq(0,1,.1),reference=1)
}
\author{
Thomas A. Gerds \email{tag@biostat.ku.dk}
}
\references{
E. Graf et al.  (1999), Assessment and comparison of
prognostic classification schemes for survival data.
Statistics in Medicine, vol 18, pp= 2529--2545.

Gerds TA, Cai T & Schumacher M (2008) The performance of
risk prediction models Biometrical Journal, 50(4), 457--479
}
\seealso{
\code{\link{pec}}
}
\keyword{survival}

