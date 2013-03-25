\name{selectFGR}
\alias{selectFGR}
\title{
Stepwise variable selection in the Fine & Gray regression competing risk model 
}
\description{
  This is a wrapper function which first selects variables in the Fine & Gray 
  regression model using \code{crrstep} from the \code{crrstep} package and 
  then returns a fitted Fine & Gray regression model with the selected variables.
}
\usage{
selectFGR(object, event, data, rule = "AIC", direction="backward")
}
\arguments{
  \item{object}{A Fine & Gray model fitted with the '\code{FGR}-function', including all the variables for stepwise selection (\code{cencode} must be 0).}
  \item{event}{Integer variable that denotes type of failure for each person, also see \code{crrstep}.}
  \item{data}{Name of a data frame containing all needed variables.}
  \item{rule}{Rule to pass on to crrstep ("AIC", "BIC" or "BICcr"), also
    see \code{crrstep} \code{crrstep}}
  \item{direction}{see \code{crrstep}}
}
\examples{
\dontrun{
library(riskRegression)
library(prodlim)
library(pec)
library(cmprsk)

d <- prodlim:::SimCompRisk(100)
newvars <- as.data.frame(matrix(runif(8*100),nrow=100))
colnames(newvars) <- c('X3','X4','X5','X6','X7','X8','X9','X10')
d <- cbind(d,newvars)

f <- FGR(formula = Hist(time, cause) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, cause=1, data = d)
sf <- selectFGR(object=f, event=cause, data=d, rule="AIC", direction="backward") 
}
}

\keyword{survival}
