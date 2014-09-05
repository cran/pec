marginal <- function(object){
  UseMethod("marginal",object)
}

marginal.default <- function(object){
  ff <- object$call$formula
  dd <- eval(object$call$data)
  fff <- reformulate("1",response=ff[[2]])
  prodlim::prodlim(fff,data=dd)
}


marginal.prodlim <- function(object){
  cc <- object$call
  ff <- cc$formula
  cc$formula <- reformulate("1",response=ff[[2]])
  eval(cc)
}

marginal.formula <- function(object){
  reformulate("1",response=object[[2]])
}
