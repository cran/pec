\name{predictEventProb}
\alias{predictEventProb}
\alias{predictEventProb.CauseSpecificCox}
\alias{predictEventProb.riskRegression}
\alias{predictEventProb.FGR}
\alias{predictEventProb.prodlim}
\title{
Predicting event probabilities (cumulative incidences) in competing risk models.
}
\description{
  Function to extract event probability predictions
  from various modeling
  approaches. The most prominent one is the combination of
  cause-specific Cox regression models which can be
  fitted with the function \code{cumincCox} from the package \code{compRisk}.
  
  The function predictEventProb is a generic function that means it invokes
  specifically designed functions depending on the 'class' of the first argument.
}
\usage{
predictEventProb(object, newdata, times, cause, ...)
\method{predictEventProb}{CauseSpecificCox}(object, newdata, times, cause, ...)
\method{predictEventProb}{riskRegression}(object, newdata, times, cause, ...)
\method{predictEventProb}{FGR}(object, newdata, times, cause, ...)
\method{predictEventProb}{prodlim}(object, newdata, times, cause, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{object}{A fitted model from which to extract predicted event probabilities}
 \item{newdata}{A data frame containing predictor variable combinations for
   which to compute predicted event probabilities.
 }
   \item{times}{A vector of times in the range of the response
    variable, for which the cumulative incidences event probabilities are computed.
  }
  \item{cause}{
    Identifies the cause of interest among the competing events.  
  }
  \item{\dots}{Additional arguments that are passed on to the current method.
  }
}
\details{
See \code{\link{predictSurvProb}}.
}
\value{  A matrix with as many rows as \code{NROW(newdata)} and as many columns
  as \code{length(times)}. Each entry should be a probability and in
  rows the values should be increasing.
}
% \references{
% }
\author{Thomas A. Gerds \email{tag@biostat.ku.dk}}
\seealso{
  See \code{\link{predictSurvProb}}.
}
% \examples{
% }
\keyword{survival}
