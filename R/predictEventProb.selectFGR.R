#########################################
# Function 'predictEventProb.selectFGR' #
#########################################

#Author: Rob C.M. van Kruijsdijk
#Date: 24-02-2013

selectFGR <- function(object,event,data,rule="AIC",direction="backward"){
  if (!require(riskRegression)) stop("This function requires library riskRegression")
  if (!require(crrstep)) stop("This function requires library crrstep")
  crrstep.form <- reformulate(paste(object$call[[2]][3]),response=all.vars(as.formula(object$call[2]))[1])
  
  capture.output(crrstep.fit <- do.call("crrstep",list(formula=crrstep.form,data=substitute(data),etype=substitute(event), 
                                       failcode=object$cause, cencode=0, direction = direction, criterion = rule, 
                                       crr.object = FALSE, trace = FALSE)))
  
  if (length(crrstep.fit$coefficients)==0){
    newform <- reformulate("1",object$call[[2]][[2]])
    newfit <- prodlim(newform,data=data) #Is this correct for competing risks (does prodlim apply CIF or KM in this case)??
  }
  else{
    newform <- reformulate(paste(rownames(crrstep.fit$coefficients),collapse="+"),response=object$call[[2]][[2]])
    newfit <- FGR(newform,data=data,cause=object$cause)
    newfit$call$formula <- newform
  }
  out <- list(fit=newfit,In=rownames(crrstep.fit$coefficients))
  out$call <- match.call()
  class(out) <- "FGRselect"
  out
}

predictEventProb.selectFGR <- function(object,newdata,times,...){
  predictEventProb(object[[1]],newdata=newdata,times=times,...)
}

###########
# TESTING #
###########

## library(riskRegression)
## library(cmprsk)

#Simulation of data (adjusted from FGR documentation):
## d <- prodlim:::SimCompRisk(100)
## newvars <- as.data.frame(matrix(runif(8*100),nrow=100))
## colnames(newvars) <- c('X3','X4','X5','X6','X7','X8','X9','X10')
## d <- cbind(d,newvars)
## evaltime <- quantile(d$time)[3]

## mod.1<- FGR(formula = Hist(time, cause) ~ X1 + X2 + X3 + X4 + X5 + X6 + X7 + X8 + X9 + X10, cause=1, data = d)
## test.1 <- selectFGR(object=mod.1,data=d,event=cause)

