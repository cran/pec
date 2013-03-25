\name{coxboost}
\alias{coxboost}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Formula interface for function \code{CoxBoost} of package \code{CoxBoost}.
}
\description{
  Formula interface for function \code{CoxBoost} of package
  \code{CoxBoost}.
}
\usage{
coxboost(formula, data, cv=TRUE, cause, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula}{
An event-history formula for competing risks of the form \code{Hist(time,status)~sex+age}
where \code{status} defines competing events and right censored
data. The code for right censored can be controlled with argument
\code{cens.code}, see man page the function \code{\link{Hist}}.
  }
  \item{data}{
    A data.frame in which the variables of formula are defined.
  }
  \item{cv}{If \code{TRUE} perform cross-validation to optimize the
    parameter \code{stepno}. This calls the function \code{cv.CoxBoost}
    whose arguments are prefix controlled, that is \code{cv.K=7} sets
    the argument \code{K} of \code{cv.CoxBoost} to \code{7}.
    If \code{FALSE} use \code{stepno}.
  }
  \item{cause}{The cause of interest in competing risk models.}
  \item{\dots}{
    Arguments passed to either \code{CoxBoost} via \code{CoxBoost.arg}
    or to \code{cv.CoxBoost} via \code{cv.CoxBoost.arg}.
  }
}
\details{
  See  \code{CoxBoost}.
 }
\value{
  See  \code{CoxBoost}.
}
\references{
  See  \code{CoxBoost}.
}
  \author{Thomas Alexander Gerds \email{tag@biostat.ku.dk} }
\seealso{
  See  \code{CoxBoost}.
}
% \examples{
% }
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{survival}
