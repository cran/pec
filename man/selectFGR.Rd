\name{selectFGR}
\alias{selectFGR}
\title{Stepwise variable selection in the Fine & Gray regression competing risk
model}
\usage{
selectFGR(formula, data, cause = 1, rule = "AIC", direction = "backward",
  ...)
}
\arguments{
  \item{formula}{A formula whose left hand side is a
  \code{Hist} object -- see \code{\link{Hist}}.  The right
  hand side specifies (a linear combination of) the
  covariates. See examples below.}

  \item{data}{A data.frame in which all the variables of
  \code{formula} can be interpreted.}

  \item{cause}{The failure type of interest. Defaults to
  \code{1}.}

  \item{rule}{Rule to pass on to crrstep ("AIC", "BIC" or
  "BICcr"), also see \code{crrstep}}

  \item{direction}{see \code{crrstep}}

  \item{\dots}{Further arguments passed to \code{crrstep}.}
}
\description{
This is a wrapper function which first selects variables in
the Fine & Gray regression model using \code{crrstep} from
the \code{crrstep} package and then returns a fitted Fine &
Gray regression model with the selected variables.
}
\examples{
\dontrun{
library(riskRegression)
library(prodlim)
library(lava)
library(cmprsk)

m <- crModel()
m <- addvar(m,c('X1','X2','X3','X4','X5','X6','X7','X8','X9','X10'))
distribution(m,c("X2","X7","X9")) <- binomial.lvm()
regression(m,eventtime1~X1+X2+X5+X9) <- c(-1,1,0.5,0.8)
set.seed(100)
d <- sim(m,100)
## full formula
ff <- Hist(time, event) ~ X1 + X2 + X3 + X4 +X5 + X6 + X7+ X8 + X9 + X10

# Fit full model with FGR
fg <- FGR(ff,cause=1,data=d)

# Backward selection based on the AIC
sfgAIC <- selectFGR(ff, data=d, rule="AIC", direction="backward")

sfgAIC$fit # Final FGR-model with selected variables

# Risk reclassification plot at time = 4
plot(predictEventProb(fg,times=4,newdata=d),
     predictEventProb(sfgAIC,times=4,newdata=d))

# Backward selection based on the BIC, while forcing
# the last two variables (X9 and X10) in the model
sfgBIC <- selectFGR(ff, data=d, rule="BIC", direction="backward",
                  scope.min=~X9+X10)

## apparent performance
pec(list(full.model=fg,selectedAIC=sfgAIC,selectedBIC=sfgBIC),
    formula=Hist(time, event)~1,
    data=d)


## bootstrap cross-validation performance
set.seed(7)
pec(list(full.model=fg,selectedAIC=sfgAIC,selectedBIC=sfgBIC),
    formula=Hist(time, event)~1,
    data=d,
    B=5,
    splitMethod="bootcv")
}
}
\author{
Rob C.M. van Kruijsdijk
\email{R.C.M.vanKruijsdijk@umcutrecht.nl}

Thomas Alexander Gerds \email{tag@biostat.ku.dk}
}
\keyword{survival}

