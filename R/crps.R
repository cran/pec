crps <- function(object,
                 who,
                 what,
                 times,
                 ## weights=NULL,
                 start){
  stopifnot(class(object)[1] == "pec")
  
  ## find the prediction models
  if (missing(who)) who <- 1:length(object$models)
  else
    if (!is.numeric(who))
      who <- names(object$models)[match(who,names(object$models))]
  
  ## times
  object.times <- object$time
  if(missing(times)) times <- object$maxtime
  if (!(object$exact || length(object.times)>100))
    warning("Only ", length(time)," time point",ifelse(length(times)==1,"","s")," used")
  
  ##  time range
  if (missing(start))
    start <- object$start
  
  ## list of values 
  if (missing(what) || is.null(what)){
    ##     print(match(c("pred.error","apparent.error","BootB0.error","boot.632.error","boot.632plus.error"),names(object)))
    ##     found.defaults <- match(c("pred.error","AppErr","OutOfBagErr","B632Err"),names(object),nomatch=0)
    ##     found.defaults <- grep(c("Err$"),names(object),val=TRUE)
    ##     what <- names(object)[found.defaults]
    what <- grep(c("Err$"),names(object),val=TRUE)
  }
  Dint <- function(x,a,b,grid){
    if ((b-a)<0)
      0
    else
      1/(b-a) * sum(x[grid>=a & grid<b]*diff(c(grid[grid>=a & grid<b],b)))
  }
  ##   if (length(weights)>0){
  ##     stopifnot(length(weights)==length(object$time))
  ##   }
  ##   else weigths <- NULL
  out <- lapply(what,function(w){
    xmat <- object[[w]][who]
    ##     if (length(weights>0)){
    ##       print("Weighted sum")
    ##       xmat <- lapply(xmat,function(u){u*weights})
    ##     }
    ##     print(xmat)
    y <- sapply(times,function(t){
      intx <- sapply(xmat, function(x){Dint(x=x,a=start,b=t,grid=object.times)})
    })
    if (!is.null(dim(y)))
      colnames(y) <- paste("t",times,sep=".")
    y
  })
  NW <- length(what)
  NT <- length(times)
  if (NW==1)
    out <- out[[1]]
  else
    names(out) <- what
  if (NT==1){
    if(NW>1){
      out <- do.call("cbind",out)
      colnames(out) <- what
    }
  }
  out
}
