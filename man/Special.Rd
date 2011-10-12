\name{Special}
\alias{Special}
\title{
  Drawing bootstrapped cross-validation curves and
  the .632 or .632plus error of models. The prediction error for an
  optional benchmark model can be added together with bootstrapped
  cross-validation error and apparent errors.
}
\description{
  This function is invoked and controlled by \code{plot.pec}. 
}
\usage{
Special(x,y,addprederr,models,bench,benchcol,times,maxboot,bootcol,col,lty,lwd)
}
\arguments{
  \item{x}{an object of class 'pec' as returned by the
    \code{pec} function.}
  \item{y}{
    Prediction error values.
  }
  \item{addprederr}{
    Additional prediction errors. The options are bootstrap
    cross-validation errors or apparent errors.
  }
  \item{models}{
    One model also specified in \code{pec} for which the \code{predErr} in
    \code{plot.pec} is to be drawn.
  }
  \item{bench}{
    A benchmark model (also specified in \code{pec}) for which the \code{predErr} in
    \code{plot.pec} is to be drawn. 
  }
  \item{benchcol}{
    Color of the benchmark curve.
  }
  \item{times}{
    Time points at which the curves must be plotted.
  }
  \item{maxboot}{
    Maximum number of bootstrap curves to be added. Default is all.
  }
  \item{bootcol}{
    Color of the bootstrapped curves. Default is 'gray77'.
  }
  \item{col}{
    Color of the different error curves for \code{models}.
}
  \item{lty}{
    Line type of the different error curves for \code{models}.
  }
  \item{lwd}{
    Line width of the different error curves for \code{models}.
  }
}
\details{
  This function should not be called directly. The arguments can be
  specified as \code{Special.arg} in the call to \code{plot.pec}.
}
\value{
 Invisible object. 
}
\seealso{
  \code{\link{plot.pec}}
}