\name{plot.pec}
\alias{plot.pec}

\title{Plotting prediction error curves}
\description{
  Plotting prediction error curves for one or more prediction models. 
}

\usage{
\method{plot}{pec}(x,
    what,
    models,
    xlim=c(x$start,x$minmaxtime),
    ylim=c(0,0.3),
    xlab="Time",
    ylab,
    axes=TRUE,
    col,
    lty,
    lwd,
    type,
    smooth=FALSE,
    add.refline=FALSE,
    add=FALSE,
    legend=ifelse(add,FALSE,TRUE),
    special=FALSE,
     ...)
}

\arguments{
  \item{x}{ Object of class \code{pec} obtained with function \code{\link{pec}}.
  }
  \item{what}{The name of the entry in \code{x}. Defauls to \code{PredErr} 
    Other choices are \code{AppErr}, \code{BootCvErr}, \code{Boot632}, \code{Boot632plus}.
  }
  \item{models}{ Specifies models in \code{x$models} for which the prediction
    error curves are drawn. Defaults to all models.
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
  \item{col}{ Vector of colors given to the curves of \code{models} in the order
    determined by \code{models}.
  }
  \item{lty}{ Vector of lty's given to the curves of \code{models} in the order
    determined by \code{models}.
  }
  \item{lwd}{ Vector of lwd's given to the curves of \code{models} in the order
    determined by \code{models}.
  }
  \item{type}{Plotting type: either \code{"l"} or \code{"s"}, see \code{lines}. 
  }
  \item{smooth}{ Logical. If \code{TRUE} the plotting values are smoothed with the
    function \code{\link{smooth}} kind="3R".
  }
  \item{add.refline}{Logical. If \code{TRUE} a dotted horizontal line is drawn as a
    symbol for the naive rule that predicts probability .5 at all
    cutpoints (i.e. time points in survival analysis).
  }
  \item{add}{Logical. If \code{TRUE} only lines are added to an existing device}
  \item{legend}{if TRUE a legend is plotted by calling the function legend. 
    Optional arguments of the function \code{legend} can be given in the
    form \code{legend.x=val} where x is the name of the argument and
    val the desired value. See also Details.
  }
  \item{special}{ Logical. If \code{TRUE} the bootstrap curves of \code{models} are
    plotted together with \code{predErr} of \code{models} by invoking the function
    \code{Special}. Optional arguments of the
    function \code{Special} can be given in the form
    \code{special.x=val} as with legend. See also Details.
  }
  \item{\dots}{ Extra arguments that are passed to \code{\link{plot}}.
  }
}

\details{
  From version 2.0.1 on the arguments legend.text, legend.args,
  lines.type, lwd.lines, specials are obsolete and only available for backward
  compatibility. Instead arguments for the invoked functions
  \code{legend}, \code{axis}, \code{Special}
  are simply specified as \code{legend.lty=2}. The specification is not
  case sensitive, thus
  \code{Legend.lty=2} or \code{LEGEND.lty=2} will have the same effect.
  The function \code{axis} is called twice, and arguments of the form
  \code{axis1.labels}, \code{axis1.at} are used for the time axis
  whereas \code{axis2.pos}, \code{axis1.labels}, etc. are used for the
  y-axis.

  These arguments are processed via  \code{\dots} of \code{plot.pec}
  and inside by using the function \code{resolveSmartArgs}.
  Documentation of these arguments can be found 
  in the help pages of the corresponding
  functions. 
}

\value{
  The (invisible) object.
}
\author{Ulla B. Mogensen  \email{ulmo@biostat.ku.dk}, Thomas A. Gerds \email{tag@biostat.ku.dk} }

\seealso{\code{\link{pec}}{\code{\link{summary.pec}}}{\code{\link{Special}}}{\code{\link{prodlim}}}}
\examples{

# simulate an artificial data frame
# with survival response and two predictors
\dontrun{
set.seed(280180)
N <- 300
X1 <- rnorm(N,mean=10,sd=5)
X2 <- rbinom(N,1,.4)
linPred <- .00001+0.2*X1+2.3*X2
T <- sapply(linPred,function(lp){rexp(n=1,exp(-lp))})
C <- rexp(n=300,exp(-mean(linPred)))
dat <- data.frame(time=pmin(T,C),status=as.numeric(T<=C),X1=X1,X2=X2)

# fit some candidate Cox models and compute the Kaplan-Meier estimate 

Models <- list("Kaplan.Meier"=survfit(Surv(time,status)~1,data=dat),
               "Cox.X1"=coxph(Surv(time,status)~X1,data=dat),
               "Cox.X2"=coxph(Surv(time,status)~X2,data=dat),
               "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=dat))


# compute the .632+ estimate of the generalization error 
set.seed(17100)
PredError.632plus <- pec(object=Models,
                         formula=Surv(time,status)~X1+X2,
                         data=dat,
                         exact=TRUE,
                         cens.model="marginal",
                         splitMethod="boot632plus",
                         B=10,
                         keep.matrix=TRUE,
                         verbose=TRUE)

# plot the .632+ estimates of the generalization error 
plot(PredError.632plus,xlim=c(0,30))

# plot the bootstrapped curves, .632+ estimates of the generalization error
# and Apparent error for the Cox model 'Cox.X1' with the 'Cox.X2' model
# as benchmark
plot(PredError.632plus,xlim=c(0,30), models="Cox.X1", special=TRUE, special.bench="Cox.X2", special.benchcol=2, special.addprederr="AppErr")
}

}

\keyword{survival}
