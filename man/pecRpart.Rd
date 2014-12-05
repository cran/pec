\name{pecRpart}
\alias{pecRpart}
\title{Predict survival based on rpart tree object}
\usage{
pecRpart(formula, data, ...)
}
\arguments{
  \item{formula}{passed to rpart}

  \item{data}{passed to rpart}

  \item{...}{passed to rpart}
}
\value{
list with three elements: ctree and call
}
\description{
Combines the rpart result with a stratified Kaplan-Meier
(prodlim) to predict survival
}
\examples{
library(prodlim)
library(rpart)
library(survival)
set.seed(50)
d <- SimSurv(50)
nd <- data.frame(X1=c(0,1,0),X2=c(-1,0,1))
f <- pecRpart(Surv(time,status)~X1+X2,data=d)
predictSurvProb(f,newdata=nd,times=c(3,8))
}

