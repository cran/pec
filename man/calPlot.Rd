\name{calPlot}
\alias{calPlot}
\title{Calibration plots for right censored data}
\usage{
calPlot(object, time, formula, data, splitMethod = "none", B = 1, M,
  outcome = c("pseudo", "prodlim"), showPseudo, pseudo.col = NULL,
  pseudo.pch = NULL, method = "nne", round = TRUE, bandwidth = NULL,
  q = 10, jack.density = 55, add = FALSE, diag = !add, legend = !add,
  axes = !add, xlim, ylim, xlab = "Predicted event probability",
  ylab = "Pseudo-observed event status", col, lwd, lty, pch, cause = 1,
  percent = TRUE, giveToModel = NULL, na.action = na.fail, cores = 1,
  verbose = FALSE, ...)
}
\arguments{
  \item{object}{A named list of prediction models, where
  allowed entries are (1) R-objects for which a
  \link{predictSurvProb} method exists (see details), (2) a
  \code{call} that evaluates to such an R-object (see
  examples), (3) a matrix with predicted probabilities
  having as many rows as \code{data} and as many columns as
  \code{times}. For cross-validation all objects in this
  list must include their \code{call}.}

  \item{time}{The evaluation time point at predicted event
  probabilities are plotted against pseudo-observed event
  status.}

  \item{formula}{A survival or event history formula. The
  left hand side is used to compute the expected event
  status. If \code{formula} is \code{missing}, try to
  extract a formula from the first element in object.}

  \item{data}{A data frame in which to validate the
  prediction models and to fit the censoring model. If
  \code{data} is missing, try to extract a data set from
  the first element in object.}

  \item{splitMethod}{Defines the internal validation
  design:

  \code{none/noPlan}: Assess the models in the give
  \code{data}, usually either in the same data where they
  are fitted, or in independent test data.

  \code{BootCv}: Bootstrap cross validation. The prediction
  models are trained on \code{B} bootstrap samples, that
  are either drawn with replacement of the same size as the
  original data or without replacement from \code{data} of
  the size \code{M}.  The models are assessed in the
  observations that are NOT in the bootstrap sample.}

  \item{B}{The number of cross-validation steps.}

  \item{M}{The size of the subsamples for
  cross-validation.}

  \item{outcome}{The method for estimating expected event
  status:

  \code{"pseudo"}: Use average pseudo-values.
  \code{"prodlim"}: Use the product-limit estimate, i.e.,
  apply the Kaplan-Meier method for right censored survival
  and the Aalen-Johansen method for right censored
  competing risks data.}

  \item{showPseudo}{If \code{TRUE} and
  \code{outcome=="pseudo"} the pseudo-values are shown as
  dots on the plot.}

  \item{pseudo.col}{Colour for pseudo-values.}

  \item{pseudo.pch}{Dot type (see par) for pseudo-values.}

  \item{method}{The method for estimating the calibration
  curve(s):

  \code{"nne"}: The expected event status is obtained in
  the nearest neighborhood around the predicted event
  probabilities.

  \code{"quantile"}: The expected event status is obtained
  in groups defined by quantiles of the predicted event
  probabilities.}

  \item{round}{If \code{TRUE} predicted probabilities are
  rounded to two digits before smoothing. This may have a
  considerable effect on computing efficiency in large data
  sets.}

  \item{bandwidth}{The bandwidth for \code{method="nne"}}

  \item{q}{The number of quantiles for
  \code{method="quantile"}.}

  \item{jack.density}{Gray scale for pseudo-observations.}

  \item{add}{If \code{TRUE} the line(s) are added to an
  existing plot.}

  \item{diag}{If \code{FALSE} no diagonal line is drawn.}

  \item{legend}{If \code{FALSE} no legend is drawn.}

  \item{axes}{If \code{FALSE} no axes are drawn.}

  \item{xlim}{Limits of x-axis.}

  \item{ylim}{Limits of y-axis.}

  \item{xlab}{Label for y-axis.}

  \item{ylab}{Label for x-axis.}

  \item{col}{Vector with colors, one for each element of
  object. Passed to \code{\link{lines}}.}

  \item{lwd}{Vector with line widths, one for each element
  of object. Passed to \code{\link{lines}}.}

  \item{lty}{lwd Vector with line style, one for each
  element of object.  Passed to \code{\link{lines}}.}

  \item{pch}{Passed to \code{\link{points}}.}

  \item{cause}{For competing risks models, the cause of
  failure or event of interest}

  \item{percent}{If TRUE axes labels are multiplied by 100
  and thus interpretable on a percent scale.}

  \item{giveToModel}{List of with exactly one entry for
  each entry in \code{object}. Each entry names parts of
  the value of the fitted models that should be extracted
  and added to the value.}

  \item{na.action}{Passed to \code{\link{model.frame}}}

  \item{cores}{Number of cores for parallel computing.
  Passed as value of argument \code{mc.cores} to
  \code{\link{mclapply}}.}

  \item{verbose}{if \code{TRUE} report details of the
  progress, e.g. count the steps in cross-validation.}

  \item{...}{Used to control the subroutines: plot, axis,
  lines, legend. See \code{\link{SmartControl}}.}
}
\value{
list with elements: time, pseudoFrame and bandwidth (NULL
for method quantile).
}
\description{
Calibration plots for risk prediction models in right
censored survival and competing risks data
}
\details{
For method "nne" the optimal bandwidth with respect to is
obtained with the function \code{\link{dpik}} from the
package \code{KernSmooth} for a box kernel function.
}
\examples{
library(prodlim)
library(lava)
library(riskRegression)
library(survival)
set.seed(13)
m <- crModel()
regression(m, from = "X1", to = "eventtime1") <- 1
regression(m, from = "X2", to = "eventtime1") <- 1
m <- addvar(m,c("X3","X4","X5"))
distribution(m, "X1") <- binomial.lvm()
distribution(m, "X4") <- binomial.lvm()
d1 <- sim(m,100)
d2 <- sim(m,100)
csc <- CSC(Hist(time,event)~X1+X2+X3+X4+X5,data=d1)
fgr <- FGR(Hist(time,event)~X1+X2+X3+X4+X5,data=d1,cause=1)
predict.crr <- cmprsk:::predict.crr
par(mar=c(5,5,5,5),cex=1.3)
calPlot(list(csc,fgr),
        time=5,
        legend.x=-0.3,
        legend.y=1.35,
        ylab="Observed event status",
        legend.legend=c("Cause-specific Cox regression","Fine-Gray regression"),
        legend.xpd=NA)
}
\author{
Thomas Alexander Gerds \email{tag@biostat.ku.dk}
}
\keyword{survival}

