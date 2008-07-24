\name{pec}
\alias{pec}
\alias{pec.list}
\alias{pec.aalen}
\alias{pec.prodlim}
\alias{pec.survfit}
\alias{pec.cph}
\alias{pec.coxph}
\alias{pec.glm}
\alias{pec.glm}

%- Also NEED an '\alias' for EACH other topic documented here.
\title{Prediction error curves}
\description{
  Evaluating the performance of risk prediction
  models in survival analysis. The Brier score is a weighted average
  of the squared distances between the observed survival status and
  the predicted survival probability of a model. Roughly the weights correspond
  to the probabilities of not being censored and they might depend on covariates.
  Prediction error curves are obtained when the Brier score is followed
  over time. Bootstrap-crossvalidation techniques are implemented to
  estimate the generalization error. 
}

\usage{
pec(object,...)
\method{pec}{list}(object,
    formula,
    data,
    times,
    start,
    maxtime,
    exact=TRUE,
    exactness=100,
    fillChar=NA,
    cens.model="cox",
    replan="none",
    B,
    M,
    model.args=NULL,
    model.parms=NULL,
    keep.matrix=FALSE,
    import=NULL,
    export=NULL,
    na.accept=0,
    verbose=TRUE,
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
    A named list of prediction models. Each entry
    is either an R-object for which a \link{predictSurvProb}
    method exists (see details) or a \code{call} that can be evaluated to such an R-object.
    For resampling (e.g. when \code{replan} is "boot632plus")
    all the models in this list must include their \code{call} in their
    value.
  }
  \item{formula}{
    A formula. If missing the formula of the
    first model in list is used.
    The left hand side is used to finde the status response
    variable in \code{data}. For right censored data, the right hand
    side of the formula is
     used to specify conditional censoring models. For
    example, set \code{Surv(time,status)~x+y} and
    \code{cens.model="cox"}. Then the weights are based on a Cox regression model for the
    censoring times with
    predictors `x' and `y'.
    Note that the usual coding is assumed:
    \code{status=0} for censored times and that each variable name that
    appears in \code{formula}
    must be the column name in \code{data}. If there are no
    covariates, i.e. \code{formula=Surv(time,status)~1} the
    \code{cens.model} is coerced to \code{"marginal"} and the
    Kaplan-Meier estimator for the censoring times is used to calculate the
    weights.
  }
  
  \item{data}{A data frame in which to validate the prediction models
    and to fit the censoring model.
    If missing the data of the
    first model in list is used.
  }
  \item{times}{A vector of timepoints. At each timepoint the prediction
    error curves
    are estimated. If missing and \code{exact==TRUE} use all the
    unique values of the response variable.
    If missing and \code{exact==FALSE} use a 
    equidistant grid of values between \code{start} and \code{maxtime}.
    The distance is determined by \code{exactness}.
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
    censoring times. Requires the timereg package maintained by Thomas Scheike. 
  }
  \item{replan}{Method for estimating the prediction error curves.
    
    \code{none}:{Assess the models in the same data where they are
      fitted. 
    }

    \code{boot}: {
      Bootstrap the prediction error curves \code{B} times WITHOUT cross-validation.
    }
    \code{cvK}: {K-fold cross-validation, i.e. \code{cv10}
      for 10-fold cross-validation.
      After splitting the data in K subsets, the prediction models (ie
      those specified in \code{object}) are evalutated on the data omitting the Kth
      subset (training step). The prediction error is estimated with the Kth subset
      (validation step).
      
      The random splitting is repeated \code{B} times and the estimated prediction error
      curves are obtained by averaging.
    }
    
    \code{OutOfBag}: {Bootstrap cross validation. The prediction models are
      trained on \code{B} bootstrap samples, that are either drawn with or
      without replacement from \code{data} of the size \code{M}.
      The models are assessed in the observations that are NOT in the
      bootstrap sample.}

    \code{Boot632}: {Linear combination of AppErr and OutOfBag using
      the constant weight .632.}

    \code{Boot632plus}: {Linear combination of AppErr and OutOfBag using
      weights dependent on how the models perform in permuted data.
    }

    \code{NoInf}: {
      Assess the models in permuted data.
    }
    
  }
  \item{B}{Number of bootstrap samples. The default depends on argument \code{replan}.
    When \code{replan in c("OutOfBag","Boot632","Boot632plus"} the
    default is 100. For \code{replan="cvK"} \code{B} is the number of cross-validation cycles,
    and -- default is 1.
    For \code{replan="none"} \code{B} is the number of bootstrap
    simulations e.g. to obtain bootstrap confidence limits -- default is 0.}
  \item{M}{The size of the bootstrap samples for resampling without replacement.
  }
  
  \item{model.args}{Experimental. List of extra arguments that
    can be passed to the \code{predictSurvProb} methods. The list must have
    an entry for each entry in \code{object}.
  }

  \item{model.parms}{Experimental.
    List of with exactly one entry for each entry in \code{object}.
    Each entry names parts of the value of the fitted models that should
    be extracted and added to the value. 
  }

  \item{keep.matrix}{
    Logical. If \code{TRUE} add all \code{B} prediction error curves to
    the output.
  }
  \item{import}{Experimental.
    Import the bootstrap matrix from an external file.
    \code{import} is a list of arguments that
    passed to internal function pecImportIndex
  }
  \item{export}{Experimental.
    Export the  bootstrap matrix to an external file.
    \code{export} is a list of arguments that
    passed to internal function pecExportIndex together
    with the matrix.
  }
  \item{na.accept}{
    Only for "OutOfBag" error.
    The maximal number of failures during training
    the models to the bootstrap samples. Usually a
    small number relative to \code{B}.
  }
  \item{verbose}{if \code{TRUE} the procedure is reporting details of
    the progress, e.g. it prints the current step in resampling
    procedures.
  }
  \item{...}{Difficult to explain}
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
    \item{}{\code{cph}, \code{psm} from \code{library(Design)}}
    \item{}{\code{prodlim} from \code{library(prodlim)}}
    \item{}{\code{glm} from \code{library(stats)}}
  }
}

\value{
  A \code{pec} object for which \code{print}, \code{summary},
  and \code{plot} methods defined.
  The object includes the following components:
  \item{PredErr}{
    The estimated prediction error according to the
    \code{replan}. A matrix where each column represents the
    estimated prediction error of a fit at
    the time points in time.
  }
  \item{AppErr}{
    The training error or apparent error obtained when
    the model(s) are evaluated in the same data where
    they were trained. Only if \code{replan} is one of
    "NoInf", "cvK", "OutOfBag", "Boot632" or "Boot632plus".
  }
  \item{OutOfBagErr}{
    The prediction error when
    the model(s) are trained in the bootstrap sample and evaluated in
    the data that are not in the bootstrap sample.
    Only if \code{replan} is one of
    "Boot632" or "Boot632plus". When \code{replan="OutOfBag"}
    then the \code{OutOfBagErr} is stored in the component \code{PredErr}.
  }
  \item{NoInfErr}{
    The prediction error when
    the model(s) are evaluated in the permuted data.
    Only if \code{replan} is one of
    "OutOfBag", "Boot632", or "Boot632plus".
    For \code{replan="NoInf"}
    the \code{NoInfErr} is stored in the component \code{PredErr}.
  }
  \item{weight}{
    The weight used to linear combine the \code{AppErr} and the
    \code{OutOfBagErr} 
    Only if \code{replan} is one of
    "Boot632", or "Boot632plus".
  }
  \item{overfit}{
    Estimated \code{overfit} of the model(s).
    See Efron \& Tibshirani (1997, Journal of the American Statistical Association) and Gerds \& Schumacher (2007, Biometrics).
    Only if \code{replan} is one of
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
  \item{method}{The method used for estimation of the prediction error.}
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

\author{Thomas Gerds \email{tag@biostat.ku.dk} }

% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{\code{\link{plot.pec}}, \code{\link{summary.pec}}, \code{\link{R2}}, \code{\link{crps}}}
\examples{

# simulate an artificial data frame
# with survival response and two predictors

set.seed(130971)
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

# compute the apparent prediction error 
PredError <- pec(object=Models,
                  formula=Surv(time,status)~X1+X2,
                  data=dat,
                  exact=TRUE,
                  cens.model="marginal",
                  replan="none",
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
                  replan="boot632plus",
                  B=100,
                  verbose=TRUE)

print(PredError.632plus,times=seq(5,30,5))
summary(PredError.632plus)
plot(PredError.632plus,xlim=c(0,30))

% ##---- Should be DIRECTLY executable !! ----
% ##-- ==>  Define data, use random,
% ##--	or do  help(data=index)  for the standard data sets.
 }
\keyword{survival}% at least one, from doc/KEYWORDS
