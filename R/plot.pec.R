# --------------------------------------------------------------------
# plot.pec function
# last update: 24 Sep 2010 (11:20)
# --------------------------------------------------------------------

# {{{ function arguments
plot.pec <- function(x,
                     what,
                     models,
                     xlim=c(x$start,x$minmaxtime),
                     ylim=c(0,0.3),
                     xlab="Time",
                     ylab,
                     axes=TRUE,
                     col,
                     lty,
                     lwd,
                     type,
                     smooth=FALSE,
                     add.refline=FALSE,
                     add=FALSE,
                     legend=ifelse(add,FALSE,TRUE),
                     special=FALSE,
                     ...){
  # }}}
# {{{ backward compatibility 

  allArgs <- match.call()
  allArgs

  # }}}
  # {{{ find the estimate to be plotted
  if (missing(what)){
    ## backward compatibility
    if (match("what",names(allArgs),nomatch=0)){
      what <- eval(allArgs$what)
    }
    else{
      # for backward compability
      if (match("PredErr",names(x),nomatch=0))
        what <- "PredErr"
      else{
        what <- switch(x$splitMethod$internal.name,
                          "noPlan"="AppErr",
                          paste(x$splitMethod$internal.name,"Err",sep=""))
      }
    }
  }
  if (0==(match(what,names(x),nomatch=0)))
    stop("Estimate \"",what,"\" not found in object")
  # }}}
  # {{{ which models
  # backward compability
  if (missing(models))
    if (match("who",names(allArgs),nomatch=0))
      models <- eval(allArgs$who)
    else
      models <- 1:length(x$models)
  if(!is.numeric(models))
    models <- names(x$model)[match(models,names(x$models))]
  
  a  <- x$time >= xlim[1]
  b <- x$time <= xlim[2]
  at <- (a & b)
  X <- x$time[at]
  y <- do.call("cbind",x[[what]][models])[at,,drop=FALSE]
   
   if (length(y)==0) stop("No plotting values: check if x[[what]][models] is a list of numeric vectors.")
    uyps <- unlist(y)
    uyps <- uyps[!is.infinite(uyps)]
    max.y <- max(uyps,na.rm=T)
    ymax <- max(max.y,ylim[2])
    if (max.y>ylim[2])
      ylim <- if (what=="PredErr")
        c(0,ceiling(ymax*10)/10)
      else
        c(0,ceiling(max(unlist(y),na.rm=T)*10))/10
  
# }}}
# {{{ Check for missings 
  nfit <- ncol(y)
  if (missing(ylab)) ylab <- "Prediction error"
  if (missing(xlab)) xlab <- "Time"
  if (missing(col)) col <- 1:nfit
  if (missing(lty)) lty <- rep(1, nfit)
  if (missing(lwd)) lwd <- rep(2, nfit)
  if (length(col)< nfit) col <- rep(col, nfit)
  if (length(lty) < nfit) lty <- rep(lty, nfit)
  if (length(lwd) < nfit) lwd <- rep(lwd, nfit)
  if (missing(type))
    if (!x$exact || smooth) type <- "l" else type <- "s"
# }}}  
# {{{ creating arguments
axis1.DefaultArgs <- list()
axis2.DefaultArgs <- list()  
plot.DefaultArgs <- list(x=0,
                         y=0,
                         type = "n",
                         ylim = ylim,
                         xlim = xlim,
                         xlab = xlab,
                         ylab = ylab)

  special.DefaultArgs <- list(x=x,y=x[[what]],addprederr=NULL,models=models,bench=FALSE,benchcol=1,times=X,maxboot=NULL,bootcol="gray77",col=rep(1,4),lty=1:4,lwd=rep(2,4))  
  if (special)
    legend.DefaultArgs <- list(legend=NULL,lwd=NULL,col=NULL,lty=NULL,cex=1.5,bty="n",y.intersp=1,x=xlim[1],xjust=0,y=(ylim+.1*ylim)[2],yjust=1)
  else
    legend.DefaultArgs <- list(legend=if(is.numeric(models)) names(x$models)[models] else models,lwd=lwd,col=col,lty=lty,cex=1.5,bty="n",y.intersp=1,x=xlim[1],xjust=0,y=(ylim+.1*ylim)[2],yjust=1)
      
# }}}
# {{{ backward compatibility

  if (match("legend.args",names(args),nomatch=FALSE)){
    legend.DefaultArgs <- c(args[[match("legend.args",names(args),nomatch=FALSE)]],legend.DefaultArgs)
    legend.DefaultArgs <- legend.DefaultArgs[!duplicated(names(legend.DefaultArgs))]
  }
  if (match("special.args",names(args),nomatch=FALSE)){
    special.DefaultArgs <- c(args[[match("special.args",names(args),nomatch=FALSE)]],special.DefaultArgs)
    special.DefaultArgs <- special.DefaultArgs[!duplicated(names(special.DefaultArgs))]
  }
  smartA <- prodlim:::SmartControl(call=list(...),
                                   keys=c("plot","special","legend","axis1","axis2"),
                                   defaults=list("plot"=plot.DefaultArgs,
                                     "special"=special.DefaultArgs,
                                     "legend"= legend.DefaultArgs,
                                     "axis1"=axis1.DefaultArgs,
                                     "axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),
                                     "axis1"=list(side=1),
                                     "axis2"=list(side=2)),
                                   ignore.case=TRUE,
                                   ignore=c("what","who"),
                                   verbose=TRUE)
  

  # }}}
  # {{{ generating an empty plot
    if (!add) {
    do.call("plot",smartA$plot)
    if (axes){
      do.call("axis",smartA$axis1)
      do.call("axis",smartA$axis2)
    }
  }
# }}}
# {{{ adding the lines 
   
   if (special==TRUE){
     if (!(x$splitMethod$internal.name=="Boot632plus"||x$splitMethod$internal.name=="Boot632"))
       stop("Plotting method 'special' requires prediction error method 'Boot632plus' or 'Boot632'")
    if (is.null(x$call$keep.matrix))
      stop("Need keep.matrix")
     do.call("Special", smartA$special)
     }
   else{
    nlines <- ncol(y)
    nix <- lapply(1:nlines, function(s) {
      lines(x = X, y = y[,s], type = type, col = col[s], lty = lty[s], lwd = lwd[s])
    })
  }

  if (add.refline) abline(h=.25,lty=3,lwd=2,col=1)
# }}}
# {{{ legend - crappy solution to legend to the option special (but works)

if(legend==TRUE && !add && !is.null(names(x$models)[models])){
     save.xpd <- par()$xpd
     par(xpd=TRUE)

     # Not very elegant solution but works:
     if (special==TRUE){
       # nameModels
       if(is.numeric(models)) nameModels <- names(x$models)[smartA$special$models] else nameModels <- smartA$special$models
       #legend.legend:
       if (is.null(smartA$legend$legend))
         if (smartA$special$bench == FALSE)
           smartA$legend$legend <- c(paste(x$method$internal.name,"-",nameModels), paste(smartA$special$addprederr, "-", nameModels))
         else{
           if (is.numeric(smartA$special$bench))
             benchName <- names(x$models)[smartA$special$bench] else benchName <- smartA$special$bench 
           if (is.null(smartA$special$addprederr))
             smartA$legend$legend <- c(paste(x$splitMethod$internal.name,"-",c(benchName, nameModels)))
           else
             smartA$legend$legend <- c(paste(x$splitMethod$internal.name,"-",c(benchName, nameModels)), paste(smartA$special$addprederr, "-", nameModels))
         }
       # legend.col
       if (is.null(smartA$legend$col))
         if (smartA$special$bench == FALSE)
           smartA$legend$col <- smartA$special$col
         else
           smartA$legend$col <- c(smartA$special$benchcol,smartA$special$col)
       if (is.null(smartA$legend$lty))
         if (smartA$special$bench == FALSE)
           smartA$legend$lty <- smartA$special$lty
         else
           smartA$legend$lty <- c(1,smartA$special$lty)
       #legend.lwd
       if (is.null(smartA$legend$lwd))
         if (smartA$special$bench == FALSE)
           smartA$legend$lwd <- smartA$special$lwd
         else
           smartA$legend$lwd <- c(smartA$special$lwd[1],smartA$special$lwd)
       # run:
       do.call("legend",smartA$legend)
     }
     
      else
        do.call("legend",smartA$legend)
     par(xpd=save.xpd)
   }


# }}}  
# {{{ returning invisible and close out
   invisible(x)
}
# }}}

