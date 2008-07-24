R2 <- function(object,
               who,
               what,
               times,
               nullModel=1){
  
  stopifnot(class(object)[1] == "pec")
  
  message("Using model ",names(object$models)[nullModel]," as reference (null model) for R^2")

  # -------------------------find the models-------------------------
  
  if (missing(who))
    who <- (1:length(object$models))[-nullModel]
  else
    if (!is.numeric(who))
      who <- match(who,names(object$models))

  # -------------------------find the estimates-------------------------

  if (missing(what) || is.null(what)){
    found.defaults <- match(c("PredErr",
                              "AppErr",
                              "OutOfBagErr",
                              "NoInfErr"),
                            names(object),
                            nomatch=0)
    what <- names(object)[found.defaults]
  }

  # ---------------------------find the times---------------------------

  object.times <- object$time
  if(missing(times)) times <- object$maxtime
  if (!(object$exact || length(object.times)>100))
    warning("Only ", length(time)," time point",ifelse(length(times)==1,"","s")," used")

  # -------------evaluate the R2 measure at specified times-------------

  nix <- lapply(what,function(e){
    if (is.null(object[[e]])) stop("No values for computing R^2")
    ref.error <- object[[e]][[nullModel]]
    out <- data.frame(do.call("cbind",lapply(1:length(who),function(w){
      rr <- 1-object[[e]][[who[w]]]/ref.error
      rr[ref.error==0] <- 0
      rr
    })))
    names(out) <- names(object$models)[who]
    cat("R^2 based on the estimate stored in ",what,":\n\n")
    
    print(cbind(time=times,RR=rbind(0,out)[1+sindex(object.times,times),,drop=FALSE]))
  })

  invisible(nix)
}
