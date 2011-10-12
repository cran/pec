\name{pec}
\alias{pec}
\alias{pec.list}
\alias{pec.aalen}
\alias{pec.prodlim}
\alias{pec.survfit}
\alias{pec.cph}
\alias{pec.coxph}
\alias{pec.glm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction error curves}
\description{
  Evaluating the performance of risk prediction
  models in survival analysis. The Brier score is a weighted average
  of the squared distances between the observed survival status and
  the predicted survival probability of a model. Roughly the weights correspond
  to the probabilities of not being censored. The weights can be
  estimated depend on covariates. Prediction error curves are
  obtained when the Brier score is followed over time.
  Cross-validation based on bootstrap resampling or bootstrap
  subsampling can be applied to assess and compare the
  predictive power of various regression modelling strategies on the
  same set of data.
}
\usage{
pec(object,...)
\method{pec}{list}(object,
                     formula,
                     data,
                     times,
                     cause,
                     start,
                     maxtime,
                     exact=TRUE,
                     exactness=100,
                     fillChar=NA,
                     cens.model="cox",
                     ipcw.refit=FALSE,
                     splitMethod="none",
                     B,
                     M,
                     reference=TRUE,
                     model.args=NULL,
                     model.parms=NULL,
                     keep.index=FALSE,
                     keep.matrix=FALSE,
                     keep.models="Call",
                     keep.residuals=FALSE,
                     keep.pvalues=FALSE,
                     noinf.permute=FALSE,
                     multiSplitTest=FALSE,
                     testIBS,
                     testTimes,
                     confInt=FALSE,
                     confLevel=0.95,
                     verbose=TRUE,
                     savePath=NULL,
                     ...)

\method{pec}{coxph}(object,...)
\method{pec}{glm}(object,...)
\method{pec}{cph}(object,...)
\method{pec}{prodlim}(object,...)
\method{pec}{survfit}(object,...)
\method{pec}{aalen}(object,...)
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
    censoring times with
    predictors x1 and x2.
    Note that the usual coding is assumed:
    \code{status=0} for censored times and that each variable name that
    appears in \code{formula}
    must be the column name in \code{data}. If there are no
    covariates, i.e. \code{formula=Surv(time,status)~1} the
    \code{cens.model} is coerced to \code{"marginal"} and the
    Kaplan-Meier estimator for the censoring times is used to calculate the
    weights.
    If \code{formula} is \code{missing}, the formula of the
    first model in list is used.
  }
  \item{data}{A data frame in which to validate the prediction models
    and to fit the censoring model.
    If \code{data} is missing, the data of the
    first model in list is used.
  }
  \item{times}{A vector of timepoints. At each timepoint the prediction
    error curves
    are estimated. If \code{exact==TRUE} the \code{times} are merged
    with all the unique values of the response variable. 
    If \code{times} is missing and \code{exact==TRUE} all the
    unique values of the response variable are used.
    If missing and \code{exact==FALSE} use a 
    equidistant grid of values between \code{start} and \code{maxtime}.
    The distance is determined by \code{exactness}.
  }
  \item{cause}{
    For competing risks, the event of interest. Defaults to the
    first state of the response, which is obtained by evaluating the
    left hand side of \code{formula} in \code{data}.
  }
  \item{start}{Minimal time for estimating the prediction error curves.
    If missing and \code{formula} defines a \code{Surv} or \code{Hist}
    object then \code{start} defaults to \code{0}, 
    otherwise to the smallest observed value of the response variable. \code{start}
    is ignored if \code{times} are given.
  }
  \item{maxtime}{Maximal time for estimating the prediction error curves. If missing
    the largest value of the response variable is used. 
  }
  \item{exact}{Logical. If \code{TRUE}
    estimate the prediction error curves at all the  unique values of
    the response variable. If \code{times} are given and
    \code{exact=TRUE} then the \code{times} are merged with the unique values of
    the response variable.
  }
  \item{exactness}{An integer that determines how many equidistant
    gridpoints are used between \code{start} and \code{maxtime}.
    The default is 100.
  }
  \item{fillChar}{Symbol used to fill-in places where the
     values of the prediction error curves are not available. The default is
    \code{NA}.}
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
    censoring times. Requires the \code{timereg} package.
  }
  \item{ipcw.refit}{
    If \code{TRUE} the inverse probability of
 censoring weigths are estimated separately in each training set during cross-validation.
  }
  \item{splitMethod}{SplitMethod for estimating the prediction error curves.
    
    \code{none/noPlan}: Assess the models in the same data where they are
      fitted. 
      \code{boot}: DEPRECIATED.
      
    \code{cvK}: K-fold cross-validation, i.e. \code{cv10}
      for 10-fold cross-validation.
      After splitting the data in K subsets, the prediction models (ie
      those specified in \code{object}) are evalutated on the data omitting the Kth
      subset (training step). The prediction error is estimated with the Kth subset
      (validation step).
      
      The random splitting is repeated \code{B} times and the estimated prediction error
      curves are obtained by averaging.
    
    \code{BootCv}: Bootstrap cross validation. The prediction models are
      trained on \code{B} bootstrap samples, that are either drawn with
      replacement of the same size as the original data or
      without replacement from \code{data} of the size \code{M}.
      The models are assessed in the observations that are NOT in the
      bootstrap sample.

    \code{Boot632}: Linear combination of AppErr and BootCvErr using
      the constant weight .632.
      
    \code{Boot632plus}: Linear combination of AppErr and BootCv using
    weights dependent on how the models perform in permuted data.

    \code{loocv}: Leave one out cross-validation.

    \code{NoInf}: Assess the models in permuted data.
  }
  \item{B}{Number of bootstrap samples. The default depends on argument \code{splitMethod}.
    When \code{splitMethod} in c("BootCv","Boot632","Boot632plus") the
    default is 100. For \code{splitMethod="cvK"} \code{B} is the number of cross-validation cycles,
    and -- default is 1.
    For \code{splitMethod="none"} \code{B} is the number of bootstrap
    simulations e.g. to obtain bootstrap confidence limits -- default is 0.}
  \item{M}{The size of the bootstrap samples for resampling without
    replacement. Ignored for resampling with replacement.
  }
  \item{reference}{Logical. If \code{TRUE} add the marginal Kaplan-Meier
  prediction model as a reference to the list of models.}

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
    index from the output list which otherwise is included in the splitMethod
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
    Logical. If \code{TRUE} keep the patient individual residuals at \code{testTimes}.
  }
  \item{keep.pvalues}{
    For \code{multiSplitTest}. If \code{TRUE} keep the pvalues from the single splits.
  }
  \item{noinf.permute}{
    If \code{TRUE} the noinformation error is approximated using permutation.
  }
  \item{multiSplitTest}{
    If \code{TRUE} the test proposed by van de Wiel et al. (2009) is
    applied. Requires subsampling bootrap cross-validation, i.e.
    that \code{splitMethod} equals \code{bootcv} and that
    \code{M} is specified.
  }
  \item{testIBS}{
    A range of time points for testing differences between models
    in the integrated Brier scores.
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
  \item{...}{Not used.}
}
\details{
  Missing data in the response or in the input matrix cause a failure. 
  
  The status of the continuous response variable at cutpoints
  (\code{times}), ie status=1 if the response value exceeds the cutpoint and
  status=0 otherwise, is compared to predicted event status probabilities
  which are provided by the prediction models on the basis of covariates.
  The comparison is done with the Brier score: the quadratic difference
  between 0-1 response status and predicted probability.
  
  There are two different sources for bias when estimating prediction
  error in right censored survival problems: censoring and high
  flexibility of the prediction model. The first is controlled by
  inverse probability of censoring weighting, the second can be
  controlled by special Monte Carlo simulation. In each step, the
  resampling procedures reevaluate the prediction model. 
  Technically this is done by replacing the argument
  \code{object$call$data}
  by the current subset or bootstrap sample of the full data.

  For each prediction model there must
  be a \code{predictSurvProb} method: for example, to assess a 
  prediction model which evaluates to a \code{myclass} object one defines 
  a function called \code{predictSurvProb.myclass} with arguments
  \code{object,newdata,cutpoints,train.data,...}

  Such a function takes the object which was fitted with train.data and
  derives a matrix with predicted event status probabilities for each subject
  in newdata (rows) and each cutpoint (column) of the response variable
  that defines an event status.

  Currently, \code{predictSurvProb} methods are available for the following R-objects:
  \describe{
    \item{}{\code{matrix}}
    \item{}{\code{aalen}, \code{cox.aalen} from \code{library(timereg)}}
    \item{}{\code{mfp} from \code{library(mfp)}}
    \item{}{\code{phnnet}, \code{survnnet} from \code{library(survnnet)}}
    \item{}{\code{rpart} (from \code{library(rpart)})}
    \item{}{\code{coxph}, \code{survfit} from \code{library(survival)}}
    \item{}{\code{cph}, \code{psm} from \code{library(rms)}}
    \item{}{\code{prodlim} from \code{library(prodlim)}}
    \item{}{\code{glm} from \code{library(stats)}}
  }
}

\value{
  A \code{pec} object. See also the help pages of the corresponding \code{print}, \code{summary},
  and \code{plot} methods.
  The object includes the following components:
  \item{PredErr}{
    The estimated prediction error according to the
    \code{splitMethod}. A matrix where each column represents the
    estimated prediction error of a fit at
    the time points in time.
  }
  \item{AppErr}{
    The training error or apparent error obtained when
    the model(s) are evaluated in the same data where
    they were trained. Only if \code{splitMethod} is one of
    "NoInf", "cvK", "BootCv", "Boot632" or "Boot632plus".
  }
  \item{BootCvErr}{
    The prediction error when
    the model(s) are trained in the bootstrap sample and evaluated in
    the data that are not in the bootstrap sample.
    Only if \code{splitMethod} is one of
    "Boot632" or "Boot632plus". When \code{splitMethod="BootCv"}
    then the \code{BootCvErr} is stored in the component \code{PredErr}.
  }
  \item{NoInfErr}{
    The prediction error when
    the model(s) are evaluated in the permuted data.
    Only if \code{splitMethod} is one of
    "BootCv", "Boot632", or "Boot632plus".
    For \code{splitMethod="NoInf"}
    the \code{NoInfErr} is stored in the component \code{PredErr}.
  }
  \item{weight}{
    The weight used to linear combine the \code{AppErr} and the
    \code{BootCvErr} 
    Only if \code{splitMethod} is one of
    "Boot632", or "Boot632plus".
  }
  \item{overfit}{
    Estimated \code{overfit} of the model(s).
    See Efron \& Tibshirani (1997, Journal of the American Statistical Association) and Gerds \& Schumacher (2007, Biometrics).
    Only if \code{splitMethod} is one of
    "Boot632", or "Boot632plus".
  }
  \item{call}{The call that produced the object} 
  \item{time}{The time points at which the prediction error curves change.}
  \item{ipcw.fit}{The fitted censoring model that was used for
    re-weigthing the Brier score residuals. See Gerds \& Schumacher
    (2006, Biometrical Journal)}
  \item{n.risk}{The number of subjects at risk for all time points.}
  \item{models}{The prediction models fitted in their own data.}
  \item{cens.model}{The censoring models.}
  \item{maxtime}{The latest timepoint where the prediction error curves are
    estimated.}
  \item{start}{The earliest timepoint where the prediction error curves are
    estimated.}
  \item{exact}{\code{TRUE} if the prediction error curves are
    estimated at all unique values of the response in the full data.}
  \item{splitMethod}{The splitMethod used for estimation of the overfitting bias.}
}

\references{E. Graf et al.  (1999),
  Assessment and comparison of prognostic classification schemes for
  survival data. Statistics in Medicine, vol 18, pp= 2529--2545.

  Efron, Tibshirani (1997) Journal of the American Statistical Association 92, 548--560
  Improvement On Cross-Validation: The .632+ Bootstrap Method.
  
  Gerds, Schumacher (2006), 
  Consistent estimation of the expected Brier score in general survival
  models with right-censored event times. Biometrical Journal, vol 48, 1029--1040.
  
  Thomas A. Gerds, Martin Schumacher (2007)
  Efron-Type Measures of Prediction Error for Survival Analysis
  Biometrics (OnlineEarly Articles).
  doi:10.1111/j.1541-0420.2007.00832.x

  Martin Schumacher, Harald Binder, and Thomas Gerds. Assessment of
  survival prediction models based on microarray data. Bioinformatics,
  23(14):1768-74, 2007.
}

\author{Thomas Alexander Gerds \email{tag@biostat.ku.dk} }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{plot.pec}}, \code{\link{summary.pec}}, \code{\link{R2}}, \code{\link{crps}}}
\examples{

# simulate an artificial data frame
# with survival response and two predictors

set.seed(130971)
dat <- SimSurv(300)

# fit some candidate Cox models and compute the Kaplan-Meier estimate 

Models <- list("Cox.X1"=coxph(Surv(time,status)~X1,data=dat,y=TRUE),
              "Cox.X2"=coxph(Surv(time,status)~X2,data=dat,y=TRUE),
              "Cox.X1.X2"=coxph(Surv(time,status)~X1+X2,data=dat,y=TRUE))

# compute the apparent prediction error 
PredError <- pec(object=Models,
                  formula=Surv(time,status)~X1+X2,
                  data=dat,
                  exact=TRUE,
                  cens.model="marginal",
                  splitMethod="none",
                  B=0,
                  verbose=TRUE)

print(PredError,times=seq(5,30,5))
summary(PredError)
plot(PredError,xlim=c(0,30))

# compute the .632+ estimate of the generalization error 
set.seed(17100)
PredError.632plus <- pec(object=Models,
                  formula=Surv(time,status)~X1+X2,
                  data=dat,
                  exact=TRUE,
                  cens.model="marginal",
                  splitMethod="Boot632plus",
                  B=100,
                  verbose=TRUE)

print(PredError.632plus,times=seq(5,30,5))
summary(PredError.632plus)
plot(PredError.632plus,xlim=c(0,30))
\dontrun{
# do the same again but now parallized
set.seed(17100)
library(doMC)
registerDoMC()
PredError.632plus <- pec(object=Models,
                  formula=Surv(time,status)~X1+X2,
                  data=dat,
                  exact=TRUE,
                  cens.model="marginal",
                  splitMethod="Boot632plus",
                  B=100,
                  verbose=TRUE)
}

# ---------------- competing risks -----------------

\dontrun{
library(survival)
library(compRisk)
library(cmprsk)
library(pec)
data(pbc)
f1  <- CRR(Hist(time,status)~sex+edema,cause=2,data=pbc)
f2  <- CRR(Hist(time,status)~sex,cause=2,data=pbc)
f3  <- compRisk(Hist(time,status)~sex+edema,cause=2,data=pbc)
f4  <- compRisk(Hist(time,status)~sex+edema,cause=2,data=pbc,link="prop")
p1 <- pec(list(f1,f2,f3,f4),formula=Hist(time,status)~1,data=pbc,cause=2)
}
}
\keyword{survival}% at least one, from doc/KEYWORDS
