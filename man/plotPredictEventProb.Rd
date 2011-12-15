\name{plotPredictEventProb}
\alias{plotPredictEventProb}
\title{ Plotting predicted survival curves.}
\description{
  Ploting time-dependent event risk predictions.
}
\usage{
plotPredictEventProb(x,
          newdata,
          times,
          cause=1,
          xlim,
          ylim,
          xlab,
          ylab,
          axes=TRUE,
          col,
          density,
          lty,
          lwd,
          add=FALSE,
          legend=TRUE,
          percent=FALSE,
           ...)
}

\arguments{
  \item{x}{ Object specifying an event risk prediction model.
  }
  \item{newdata}{ A data frame with the same variable names as those that
    were used to fit the model \code{x}.
  }
  \item{times}{ Vector of times at which to return the estimated probabilities.
  }
  \item{cause}{Show predicted risk of
    events of this cause
  }
  \item{xlim}{ Plotting range on the x-axis.
  }
  \item{ylim}{ Plotting range on the y-axis.
  }
  \item{xlab}{ Label given to the x-axis.
  }
  \item{ylab}{ Label given to the y-axis.
  }
  \item{axes}{ Logical. If \code{FALSE} no axes are drawn.
  }
  \item{col}{ Vector of colors given to the survival curve.
  }
  \item{density}{Densitiy of the color -- useful for showing many (overlapping) curves.
  }
  \item{lty}{ Vector of lty's given to the survival curve.
  }
  \item{lwd}{ Vector of lwd's given to the survival curve.
  }
  \item{add}{ Logical. If \code{TRUE} only lines are added to an existing device}
  \item{legend}{ Logical. If TRUE a legend is plotted by calling the function legend. 
    Optional arguments of the function \code{legend} can be given in the
    form \code{legend.x=val} where x is the name of the argument and
    val the desired value. See also Details.
  }
  \item{percent}{Logical. If \code{TRUE} the y-axis is labeled in percent.
  }
  \item{\dots}{Parameters that are filtered by \code{\link{SmartControl}} and
    then passed to the functions: \code{\link{plot}}, \code{\link{axis}},
    \code{\link{legend}}.
  }
}

\details{
  Arguments for the invoked functions \code{legend} and \code{axis}
  are simply specified as \code{legend.lty=2}. The specification is not
  case sensitive, thus \code{Legend.lty=2} or \code{LEGEND.lty=2} will
  have the same effect. The function \code{axis} is called twice, and
  arguments of the form \code{axis1.labels}, \code{axis1.at} are used
  for the time axis whereas \code{axis2.pos}, \code{axis1.labels}, etc.
  are used for the y-axis.

  These arguments are processed via  \code{\dots} of \code{plotPredictEventProb}
  and inside by using the function \code{SmartControl}.
}

\value{
  The (invisible) object.
}
\author{Ulla B. Mogensen  \email{ulmo@biostat.ku.dk}, Thomas A. Gerds \email{tag@biostat.ku.dk} }
\seealso{{\code{\link{predictEventProb}}}{\code{\link{prodlim}}}}

\examples{

# generate some competing risk data
}
\keyword{survival}
