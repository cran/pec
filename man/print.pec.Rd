\name{print.pec}
\alias{print.pec}
\alias{summary.pec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Printing a `pec' (prediction error curve) object.}
\description{
  Print the important arguments of call and the prediction error values
  at selected time points. 
}
\usage{
\method{print}{pec}(x, times, digits = 3,what, ...)
\method{summary}{pec}(object,times,what,models,digits=3,print=TRUE,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{Object of class \code{pec}}
  \item{object}{Object of class \code{pec}}
  \item{times}{Time points at which to show the values of the prediction error curve(s)}
  \item{what}{What estimate of the prediction error curve to
    show. Should be a string matching an element of the object. The
    default is determined by splitMethod.}
  \item{models}{Which models in the list \code{object$models} should be
    shown. Defaults to all models.}
  \item{digits}{Number of decimals used in tables.}
  \item{print}{Set to FALSE to suppress printing.}
  \item{\dots}{Not used} }
\value{
  The first argument in the invisible cloak.
}
% \references{ ~put references to the literature/web site here ~ }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk} }
% \note{ ~~further notes~~ 

%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{pec}}}
\keyword{survival}% at least one, from doc/KEYWORDS
