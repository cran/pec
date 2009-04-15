\name{plot.pec}
\alias{plot.pec}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Plotting prediction error curves}
\description{
  Ploting prediction error curves for one or more
  prediction models. 
}
\usage{
\method{plot}{pec}(x,
	 what="PredErr",
	 who,
	 xlim=c(x$start,x$maxtime),
	 ylim=c(0,0.3),
	 xlab="Time",
	 ylab,
	 lwd.lines=2,
	 axes=TRUE,
	 col,
	 lty,
	 lines.type,
	 smooth=FALSE,
	 add.refline=FALSE,
	 add=FALSE,
	 legend=ifelse(add,FALSE,TRUE),
	 legend.text,
	 legend.args=list(),
	 specials=NULL,
	 ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{object of class \code{pec} obtained with function \code{\link{pec}}}
  \item{what}{the name of the entry in \code{x}
Defauls to \code{PredErr} 
Other choices are \code{AppErr}, \code{OutOfBagErr}, \code{NoInfErr}, \code{overfit}
}
\item{who}{Specifies models in \code{x$models} for which the prediction
  error curves are drawn. Defaults to all models.
}
\item{xlim}{Plotting range on the x-axis.}
\item{ylim}{Plotting range on the y-axis}
\item{xlab}{Label given to the x-axis.}
\item{ylab}{Label given to the y-axis.}
\item{lwd.lines}{Extra parameter passed to function \code{\link{lines}}.}
\item{axes}{if FALSE no axes are drawn}
\item{col}{vector of colors given to the curves of models in the order determined by \code{who}}
\item{lty}{vector of lty's given to the curves of models in the order determined by \code{who}}
\item{lines.type}{passed to \code{\link{lines}}}
\item{smooth}{if TRUE the plotting values are smoothed with the
  function \code{\link{smooth}} kind="3R"}
\item{add.refline}{if TRUE a dotted horizontal line is drawn as a
  symbol for the naive rule that predicts probability .5 at all
  cutpoints (i.e. time points in survival analysis).}
\item{add}{TRUE means that only lines are added to an existing device}
\item{legend}{TRUE means draw the legend}
\item{legend.text}{legend text in the order determined by \code{who}}
\item{legend.args}{a list with arguments that are passed to
  \code{\link{legend}}.}
\item{specials}{Experimental}
\item{\dots}{extra arguments that are passed to \code{\link{plot}}.}
}
% \details{
% }
% \value{
%   ~Describe the value returned
%   If it is a LIST, use
%   \item{comp1 }{Description of 'comp1'}
%   \item{comp2 }{Description of 'comp2'}
%   ...
% }
% \references{ ~put references to the literature/web site here ~ }
\author{Thomas Gerds \email{tag@biostat.ku.dk} }
% \note{ ~~further notes~~ 
%  ~Make other sections like Warning with \section{Warning }{....} ~
% }
\seealso{\code{\link{pec}}}
% \examples{
% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.

\keyword{survival}
