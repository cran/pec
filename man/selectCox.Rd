\name{selectCox}
\alias{selectCox}
\title{
Backward variable selection in the Cox regression model 
}
\description{
  This is a wrapper function which first selects variables in the Cox
  regression model using \code{fastbw} from the \code{rms} package and
  then returns a fitted Cox regression model with the selected variables.
}
\usage{
selectCox(formula,data,rule="aic")
}
\arguments{
  \item{formula}{A formula object with a \code{Surv} object on the
left-hand side and all the variables on the right-hand side. 
}
  \item{data}{Name of an data frame containing all  needed variables.}
  \item{rule}{The method for selecting variables. See
    \code{\link{fastbw}} for details.}
}
\details{
This function first calls \code{cph} then \code{fastbw} and finally
\code{cph} again.
}
\references{
  Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012). Evaluating Random Forests
  for Survival Analysis Using Prediction Error Curves. Journal of Statistical Software,
  50(11), 1-23. URL http://www.jstatsoft.org/v50/i11/.
}
\examples{
\dontrun{
data(GBSG2)
f <- selectCox(Surv(time,cens)~horTh + age + menostat + tsize + tgrade + pnodes + progrec + estrec ,data=GBSG2)
}
}
\keyword{survival}
