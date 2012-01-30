R2 <- function(object,
               models,
               what,
               times,
               reference=1){
  
  stopifnot(class(object)[1] == "pec")
  
  # {{{find the prediction models
  
  if (missing(models))
    models <- (1:length(object$models))[-reference]
  else
    if (!is.numeric(models))
      models <- match(models,names(object$models))
  # }}}
  # {{{ what errors
  if (missing(what) || is.null(what)){
    what <- grep(c("Err$"),names(object),value=TRUE)
  }
  # }}}
  # {{{ find the times

  object.times <- object$time
  if(missing(times)) times <- object$maxtime
  if (!(object$exact || length(object.times)>100))
    warning("Only ", length(time)," time point",ifelse(length(times)==1,"","s")," used")
  # }}}
  # {{{ for each element of what: evaluate R2 at times

  out <- lapply(what,function(e){
    if (is.null(object[[e]])) stop("No values for computing R^2")
    ref.error <- object[[e]][[reference]]
    out <- data.frame(do.call("cbind",lapply(1:length(models),function(w){
      rr <- 1-object[[e]][[models[w]]]/ref.error
      rr[ref.error==0] <- 0
      rr
    })))
    names(out) <- names(object$models)[models]
    ## cat("R^2 based on the estimate stored in ",what,":\n\n")
    ## print(cbind(time=times,RR=rbind(0,out)[1+sindex(object.times,times),,drop=FALSE]))
    cbind(time=times,RR=rbind(0,out)[1+sindex(object.times,times),,drop=FALSE])
  })
  # }}}
  # {{{ prepare output
  NW <- length(what)
  NT <- length(times)
  names(out) <- what
  # }}}
  attr(out,"reference") <- names(object$models)[reference]
  class(out) <- "R2"
  out
}
