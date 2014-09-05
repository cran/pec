\name{resolvesplitMethod}
\alias{resolvesplitMethod}
\title{Resolve the splitMethod for estimation of prediction performance}
\usage{
resolvesplitMethod(splitMethod, B, N, M)
}
\arguments{
  \item{splitMethod}{String that determines the splitMethod
  to use. Available splitMethods are none/noPlan (no
  splitting), bootcv or outofbag (bootstrap
  cross-validation), cvK (K-fold cross-validation, e.g.
  cv10 gives 10-fold), boot632, boot632plus or boot632+,
  loocv (leave-one-out)}

  \item{B}{The number of repetitions.}

  \item{N}{The sample size}

  \item{M}{For subsampling bootstrap the size of the
  subsample. Note M<N.}
}
\value{
A list with the following components \item{name }{the
official name of the splitMethod} \item{internal.name }{the
internal name of the splitMethod} \item{index}{a matrix of
indices with B columns and either N or M rows, dependent on
splitMethod} \item{B}{the value of the argument B}
\item{N}{the value of the argument N} \item{M}{the value of
the argument M}
}
\description{
The function computes a matrix of random indices obtained
by drawing from the row numbers of a data set either with
or without replacement.  The matrix can be used to
repeatedly set up independent training and validation sets.
}
\examples{
# BootstrapCrossValidation: Sampling with replacement
  resolvesplitMethod("BootCv",N=10,B=10)

  # 10-fold cross-validation: repeated 2 times
  resolvesplitMethod("cv10",N=10,B=2)

  # leave-one-out cross-validation
  resolvesplitMethod("loocv",N=10)

  resolvesplitMethod("bootcv632plus",N=10,B=2)
}
\author{
Thomas Alexander Gerds \email{tag@biostat.ku.dk}
}
\keyword{prediction}

