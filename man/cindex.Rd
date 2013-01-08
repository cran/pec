\name{cindex}
\alias{cindex}
\alias{cindex.list}
\title{Concordance index for right censored survival time data}
\description{
In survival analysis, a pair of patients is called concordant if the
risk of the event predicted by a model is lower for the patient who
experiences the event at a later timepoint. The concordance probability
(C-index) is the frequency of concordant pairs among all pairs of
subjects. It can be used to measure and compare the discriminative power
of a risk prediction models. The function provides an inverse of the probability of
censoring weigthed estimate of the concordance probability to adjust for
right censoring. Cross-validation based on bootstrap resampling or
bootstrap subsampling can be applied to assess and compare the
discriminative power of various regression modelling strategies on the
same set of data.
}
\usage{
cindex(object,...)
\method{cindex}{list}(object,
                        formula,
                        data,
                        eval.times,
                        pred.times,
                        cause,
                        cens.model="marginal",
                        ipcw.refit=FALSE,
                        ipcw.limit,
                        tiedPredictionsIn=TRUE,
                        tiedOutcomeIn=TRUE,
                        tiedMatchIn=TRUE,
                        splitMethod="noPlan",
                        B,
                        M,
                        model.args=NULL,
                        model.parms=NULL,
                        keep.models="Call",
                        keep.residuals=FALSE,
                        keep.pvalues=FALSE,
                        keep.weights=FALSE,
                        keep.index=FALSE,
                        keep.matrix=FALSE,
                        multiSplitTest=FALSE,
                        testTimes,
                        confInt=FALSE,
                        confLevel=0.95,
                        verbose=TRUE,
                        savePath=NULL,
                        slaveseed=NULL,
                        ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
    A named list of prediction models, where allowed entries are
    (1) R-objects for which a \link{predictSurvProb}
    method exists (see details), (2) a \code{call} that evaluates
    to such an R-object (see examples), (3) a matrix with predicted
    probabilities having as many rows as \code{data} and as many columns
    as \code{times}. For cross-validation all objects in this list
    must include their \code{call}.
  }
  \item{formula}{
    A survival formula. The left hand side is used to finde the status response
    variable in \code{data}. For right censored data, the right hand
    side of the formula is
     used to specify conditional censoring models. For
     example, set \code{Surv(time,status)~x1+x2} and
    \code{cens.model="cox"}. Then the weights are based on a Cox regression model for the
    censoring times with predictors x1 and x2.
    Note that the usual coding is assumed:
    \code{status=0} for censored times and that each variable name that
    appears in \code{formula}
    must be the column name in \code{data}. If there are no
    covariates, i.e. \code{formula=Surv(time,status)~1} the
    \code{cens.model} is coerced to \code{"marginal"} and the
    Kaplan-Meier estimator for the censoring times is used to calculate the
    weights.
    If formula is \code{missing} the formula of the
    first model in list is used.
  }
  \item{data}{A data frame in which to validate the prediction models
    and to fit the censoring model.
    If \code{data} is missing, the data of the
    first model in list is used.
  }
 \item{eval.times}{A vector of timepoints for evaluating the
discriminative ability. At each timepoint the c-index
is computed using only those pairs where one of the event times is known
to be earlier than this timepoint. If \code{eval.times} is \code{missing} or
\code{Inf} then the largest uncensored event time is used. 
}
 \item{pred.times}{A vector of timepoints for evaluating the prediction
models. This should either be exactly one timepoint used for all
\code{eval.times}, or be as long as \code{eval.times}, in which case the
predicted order of risk for the jth entry of \code{eval.times}
is based on the jth entry of \code{pred.times}
corresponding 
 }
  \item{cause}{
    For competing risks, the event of interest. Defaults to the
    first state of the response, which is obtained by evaluating the
    left hand side of \code{formula} in \code{data}.
  }
  \item{cens.model}{
    Method for estimating inverse probability of
    censoring weigths:
    
    \code{cox}: A semi-parametric Cox proportional hazard model is
    fitted to the censoring times
    
    \code{marginal}: The Kaplan-Meier estimator for the censoring times
    
    \code{nonpar}: Nonparametric extension of the 
    Kaplan-Meier for the censoring times using symmetric nearest
    neighborhoods -- available for
    arbitrary many strata variables on the right hand side of
    argument \code{formula} but at most one continuous
    variable. See the documentation of the functions \code{prodlim} and \code{neighborhood} from
    the prodlim package.

    \code{aalen}: The nonparametric Aalen additive model fitted to the
    censoring times. Requires the timereg package maintained by Thomas Scheike. 
  }
  \item{ipcw.refit}{
    If \code{TRUE} the inverse probability of
    censoring weigths are estimated separately in each training set during cross-validation.
  }
  \item{ipcw.limit}{
    Value between 0 and 1 (but no equal to 0!) used to cut for small
    weights in order to stabilize the estimate at late times were few
    individuals are observed.
  }
  \item{tiedPredictionsIn}{If \code{FALSE} pairs with identical predictions
    are excluded, unless also the event times are identical and uncensored and
    \code{tiedMatchIn} is set to \code{TRUE}.
  }
  \item{tiedOutcomeIn}{
    If \code{TRUE} pairs with identical and uncensored event times are excluded, unless
    also the predictions are identical and \code{tiedMatchIn} is set to
    \code{TRUE}.
  }
  \item{tiedMatchIn}{
    If \code{TRUE} then pairs with identical predictions and identical and
    uncensored event
    times are counted as concordant pairs.
  }
  \item{splitMethod}{SplitMethod for estimating the prediction error curves.
    
    \code{none/noPlan}: Assess the models in the give \code{data},
    usually either in the same data where they are
    fitted, or in independent test data. 
    
    \code{BootCv}: Bootstrap cross validation. The prediction models are
      trained on \code{B} bootstrap samples, that are either drawn with
      replacement of the same size as the original data or
      without replacement from \code{data} of the size \code{M}.
      The models are assessed in the observations that are NOT in the
      bootstrap sample.
      
    \code{Boot632}: Linear combination of AppCindex and OutOfBagCindex using
      the constant weight .632.
      
  }
  \item{B}{Number of bootstrap samples. The default depends on argument \code{splitMethod}.
    When \code{splitMethod} in c("BootCv","Boot632") the
    default is 100.
    For \code{splitMethod="none"} \code{B} is the number of bootstrap
    simulations e.g. to obtain bootstrap confidence limits -- default is
    0.}
  \item{M}{The size of the bootstrap samples for resampling without
    replacement. Ignored for resampling with replacement.
  }

  \item{model.args}{List of extra arguments that
    can be passed to the \code{predictSurvProb} methods. The list must have
    an entry for each entry in \code{object}.
  }

  \item{model.parms}{Experimental.
    List of with exactly one entry for each entry in \code{object}.
    Each entry names parts of the value of the fitted models that should
    be extracted and added to the value. 
  }

  \item{keep.index}{
    Logical. If \code{FALSE} remove the bootstrap or cross-validation
    index from the output list which otherwise is included in the method
    part of the output list. 
  }

  \item{keep.matrix}{
    Logical. If \code{TRUE} add all \code{B} prediction error curves
    from bootstrapping or cross-validation to the output.
  }
  \item{keep.models}{
    Logical. If \code{TRUE} keep the models in object. If \code{"Call"}
    keep only the \code{call} of these models. 
  }
  \item{keep.residuals}{
Experimental.
  }
  \item{keep.pvalues}{
Experimental.
  }
  \item{keep.weights}{
Experimental.
  }
  \item{multiSplitTest}{
  Experimental.
  }
  \item{testTimes}{
    A vector of time points for testing differences between models
    in the time-point specific Brier scores.
  }
  \item{confInt}{
    Experimental.
  }
  \item{confLevel}{
    Experimental.
  }
  \item{verbose}{if \code{TRUE} report details of
    the progress, e.g. count the steps in cross-validation.
  }
  \item{savePath}{Place in your filesystem (directory) where training models fitted during
    cross-validation are saved. If \code{missing} training models are not saved.
  }
    \item{slaveseed}{
    Vector of seeds, as long as \code{B}, to be given to the slaves in parallel computing.
  }
  \item{...}{Not used.}
}
\details{
Pairs with identical observed times, where one is uncensored and one is censored, are
always considered usuable (independent of the value of
\code{tiedOutcomeIn}), as it can be assumed that the event occurs at a
later timepoint for the censored observation.

For uncensored response the result equals the one obtained with the
functions \code{rcorr.cens} and \code{rcorrcens} from the \code{Hmisc}
package (see examples).
}
\value{
Estimates of the C-index. 
}
\references{

Gerds, TA and Kattan, M and Schumacher, M and Yu, C (2010)
Estimating a time-dependent concordance index for survival prediction models with covariate dependent censoring
Research report 10/7. Department of Biostatistics, University of
Copenhagen
}
\author{Thomas A Gerds \email{tag@biostat.ku.dk}}


\examples{
\dontrun{
 # simulate data based on Weibull regression  
library(prodlim)
 set.seed(13)
 dat <- SimSurv(300)
 # fit three different Cox models and a random survival forest
 # note: low number of trees for the purpose of illustration 
 library(survival)
 library(randomSurvivalForest)
 #rsf.default=randomSurvivalForest:::rsf.default
 cox12 <- coxph(Surv(time,status)~X1+X2,data=dat)
 cox1 <- coxph(Surv(time,status)~X1,data=dat)
 cox2 <- coxph(Surv(time,status)~X2,data=dat)
 rsf1 <- rsf(Survrsf(time,status)~X1+X2,data=dat,ntree=15,forest=TRUE)
 #
 # compute the apparent estimate of the C-index at different time points
 #
ApparrentCindex  <- cindex(list("Cox X1"=cox1,
		       "Cox X2"=cox2,
		       "Cox X1+X2"=cox12,
                       "RSF"=rsf1),
		  formula=Surv(time,status)~X1+X2,
		  data=dat,
		  eval.times=seq(5,500,50))
  print(ApparrentCindex)
  plot(ApparrentCindex)
 #
 # compute the bootstrap-crossvalidation estimate of the C-index at different time points
 #
set.seed(142)
bcvCindex  <- cindex(list("Cox X1"=cox1,
		       "Cox X2"=cox2,
		       "Cox X1+X2"=cox12,
                       "RSF"=rsf1),
		  formula=Surv(time,status)~X1+X2,
		  data=dat,
                  splitMethod="bootcv",
                  B=10,
		  eval.times=seq(5,500,50))
  print(bcvCindex)
  plot(bcvCindex)
}
\dontrun{
 # for uncensored data the results are the same as obtained with the function rcorr.cens from Hmisc
library(Hmisc)
set.seed(16)
dat <- SimSurv(30,cens=FALSE)
fit12 <- coxph(Surv(time,status)~X1+X2,data=dat)
fit1 <- coxph(Surv(time,status)~X1,data=dat)
fit2 <- coxph(Surv(time,status)~X2,data=dat)
Cpec <- cindex(list("Cox X1+X2"=fit12,"Cox X1"=fit1,"Cox X2"=fit2),formula=Surv(time,status)~1,data=dat,eval.times=Inf)
p1 <- predictSurvProb(fit1,newdata=dat,times=100)
p2 <- predictSurvProb(fit2,newdata=dat,times=100)
p12 <- predictSurvProb(fit12,newdata=dat,times=100)
harrelC1 <- rcorr.cens(p1,with(dat,Surv(time,status)))
harrelC2 <- rcorr.cens(p2,with(dat,Surv(time,status)))
harrelC12 <- rcorr.cens(p12,with(dat,Surv(time,status)))
harrelC1[["C Index"]]==Cpec$AppCindex[["Cox.X1"]]
harrelC2[["C Index"]]==Cpec$AppCindex[["Cox.X2"]]
harrelC12[["C Index"]]==Cpec$AppCindex[["Cox.X1.X2"]]
 #
 # competing risks 
 #
# require(lava)
m <- lvm()
X <- paste("X",seq(5),sep="")

}
}

\keyword{survival}
