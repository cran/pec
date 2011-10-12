calPlot <- function(object,
                    predTime,
                    formula,
                    data,
                    method="snne",
                    bandwidth=NULL,
                    verbose=FALSE,
                    add=FALSE,
                    showPseudo=TRUE,
                    diag=TRUE,
                    legend=TRUE,
                    axes=TRUE,
                    background=FALSE,
                    Grid=background,
                    xlim,
                    ylim,
                    xlab,
                    ylab,
                    col,
                    lwd,
                    lty,
                    pch,
                    cause=1,
                    percent=TRUE,
                    ...){
  
  # {{{ find number of objects and lines
  cobj=class(object)
  if (cobj!="list"){
    object <- list(object)
  }
  if (is.null(names(object))){
    names(object) <- sapply(object,function(o)class(o)[1])
  }
  else{
    names(object)[(names(object)=="")] <- sapply(object[(names(object)=="")],function(o)class(o)[1])
  }
  names(object) <- make.names(names(object),unique=TRUE)
  nlines <- length(object)
  # }}}
  # {{{ lines types
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- 1:nlines
  if (missing(lty)) lty <- rep(1, nlines)
  if (missing(pch)) pch <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  if (length(pch) < nlines) pch <- rep(pch, nlines)
  # }}}
  # {{{ data & formula
  if (missing(data)){
    data <- eval(object[[1]]$call$data)
    if (match("data.frame",class(data),nomatch=0)==0)
      stop("Argument data is missing.")
    else
      if (verbose)
        warning("Argument data is missing. I use the data from the call to the first model instead.")
  }
  if (missing(formula)){
    formula <- eval(object[[1]]$call$formula)
    if (match("formula",class(formula),nomatch=0)==0)
      stop("Argument formula is missing.")
    else if (verbose)
      warning("Argument formula is missing. I use the formula from the call to the first model instead.")
  }
  m <- model.frame(formula,data,na.action=na.fail)
  response <- model.response(m)
  if (match("Surv",class(response),nomatch=FALSE))
    model.type <- "survival"
  else
    model.type <- attr(response,"model")
  neworder <- order(response[,"time"],-response[,"status"])
  response <- response[neworder,,drop=FALSE]
  Y <- response[,"time"]
  status <- response[,"status"]
  data <- data[neworder,]
  # }}}
  # {{{ prediction timepoint 
  if (missing(predTime))
    predTime <- median(Y)
  else
    if (length(predTime)>1)
      stop("Please specify only one time point predTime.")
  # }}}
  # {{{ compute pseudo-values
  #  require(pseudo)
  #  jack=pseudosurv(time=Y,event=status,tmax=predTime)[[3]]
  if (model.type=="competing.risks"){
    predictHandlerFun <- "predictEventProb"
  }
  else{
    predictHandlerFun <- "predictSurvProb"
  }
  margForm <- reformulate("1",response=formula[[2]])
  margFit <- prodlim(margForm,data=data)
  jack <- jackknife(margFit,times=predTime)
  # }}}
  # {{{ call smartControls
  axis1.DefaultArgs <- list(side=1,las=1)
  axis2.DefaultArgs <- list(side=2,las=2)
  legend.DefaultArgs <- list(legend=names(object),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="topleft")
  Grid.DefaultArgs <- list(vertical=NULL,horizontal=NULL,col=c("gray88","gray99"))
  background.DefaultArgs <- list(col="gray77")
  lines.DefaultArgs <- list(type="l")
  if (missing(ylim)){
    if (showPseudo)
      ylim <- c(min(jack),max(jack))
    else
      ylim <- c(0,1)
  }
  if (missing(xlim)){
    xlim <- c(0,1)
  }
  plot.DefaultArgs <- list(x=0,y=0,type = "n",ylim = ylim,xlim = xlim,xlab = "Predicted probabilities",ylab = "Observed probabilities")
  smartA <- prodlim:::SmartControl(call= list(...),
                                   keys=c("plot","lines","legend","Grid","background","axis1","axis2"),
                                   ignore=NULL,
                                   ignore.case=TRUE,
                                   defaults=list("plot"=plot.DefaultArgs,"lines"=lines.DefaultArgs,"legend"=legend.DefaultArgs,"Grid"=Grid.DefaultArgs,"background"=background.DefaultArgs,"axis1"=axis1.DefaultArgs,"axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),"axis1"=list(side=1)),
                                   verbose=TRUE)
  # }}}
  # {{{ plot an empty frame
  if (add==FALSE){
    do.call("plot",smartA$plot)
    if (axes){
      if (percent){
        if (is.null(smartA$axis2$at)){
          smartA$axis2$at <- seq(0,1,.25)
          smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
        }
        else{
          smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
        }
        if (is.null(smartA$axis1$at)){
          smartA$axis1$at <- seq(0,1,.25)
          smartA$axis1$labels <- paste(100*smartA$axis1$at,"%")
        }
        else{
          smartA$axis1$labels <- paste(100*smartA$axis1$at,"%")
        }
      }
      do.call("axis",smartA$axis2)
      do.call("axis",smartA$axis1)
      do.call("axis",smartA$axis2)
    }
    if (background)
      do.call("prodlim:::background",smartA$background)
    if (Grid){
      if (is.null(smartA$Grid$horizontal) && !is.null(smartA$axis2$at))
        smartA$Grid$horizontal <- smartA$axis2$at
      if (is.null(smartA$Grid$vertical) && !is.null(smartA$axis1$at))
        smartA$Grid$vertical <- smartA$axis1$at
      do.call("prodlim:::Grid",smartA$Grid)
    }
  }
  if (diag)
    segments(x0=0,y0=0,x1=1,y1=1,col="gray77",lwd=2,xpd=FALSE)
  ##   do.call("abline",c(list(a=0,b=1),list(col="gray77",lwd=2,xpd=FALSE)))
  
  # }}}
  # {{{ add one smoothed calibration line for each model
  method <- match.arg(method,c("fixed","snne","loess","kernel"))
  predictions <- lapply(1:nlines,function(f){
    if (model.type=="competing.risks"){
      p=do.call(predictHandlerFun,list(object[[f]],newdata=data,times=predTime,cause=cause))
    }
    else{
      p=do.call(predictHandlerFun,list(object[[f]],newdata=data,times=predTime))
    }
    if (showPseudo) points(p,jack,col="gray")
    switch(method,
           "fixed"={
             ## x <- p
             ## y <- sapply(unique(p),function(prob){
             ## nbh <- (1:NROW(p)) prob
             ## jack
             ## }
             ## browser()
           },
           "snne"={
             plotFrame=meanNeighbors(x=p,y=jack,bandwidth=bandwidth)
             lines(plotFrame$uniqueX,plotFrame$averageY,col=col[f],lty=lty[f],lwd=lwd[f])
           },
           "loess"={
             loess.args=NULL
             loess.DefaultArgs <- list(family="symmetric",control=loess.control(iterations=0),span=.3)
             loess.args <- c(formula=jack~pp,loess.args,loess.DefaultArgs)
             loess.args <- loess.args[!duplicated(names(loess.args))]
             pp=sort(p)
             smoothJack=do.call("loess",loess.args)
             plotFrame=data.frame(x=pp,y=smoothJack$fitted)
             with(na.omit(plotFrame),lines(x,y,col=col[f],lwd=lwd[f],lty=lty[f]))
           },
           "kernel"={
             ##              require(kknn)
             message("kernel smoothing not yet implemented")
           })
    # lines(pseudoFrame$uniqueX,smooth(pseudoFrame$averageY),col=col[f],lwd=lwd)
    c(p)
  })
  # }}}
  # {{{ legend
  ## if (missing(legend)) legend=ifelse(length(object)==1,FALSE,TRUE)
  ## if (missing(legend.legend)) legend.legend=names(object)
  do.call("legend",smartA$legend)
  ## if (legend)
  ## legend(0,1,legend=legend.legend,lwd=lwd,col=col,bty="n")
  # }}}
  # {{{ invisibly output the jackknife pseudo-values
  pseudoFrame=data.frame(Y=jack)
  pseudoFrame=cbind(pseudoFrame,do.call("cbind",predictions))
  names(pseudoFrame) <- c("PseudoResponse",paste("predictProbModel",1:nlines,sep="."))
  out <- list(predTime=predTime,pseudoFrame=pseudoFrame)
  invisible(out)
  # }}}
}
