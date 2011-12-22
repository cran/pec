plotPredictEventProb <- function(x,
                                 newdata,
                                 times,
                                 cause=1,
                                 xlim,
                                 ylim,
                                 xlab,
                                 ylab,
                                 axes=TRUE,
                                 col,
                                 density,
                                 lty,
                                 lwd,
                                 add=FALSE,
                                 legend=TRUE,
                                 percent=FALSE,
                                 ...){
  # {{{ call argument

  allArgs <- match.call()
  
  # }}}
  # {{{ find times

  if(missing(times)){
    # formula
    formula <- eval(x$call$formula)
    if (match("formula",class(formula),nomatch=0)==0)
      stop("Argument formula is missing.")
    # find data
    data <- eval(x$call$data)
    # extract response
    m <- model.frame(formula,data,na.action=na.fail)
    response <- model.response(m)
    # ordering time 
    neworder <- order(response[,"time"],-response[,"status"])
    response <- response[neworder,,drop=FALSE]
    times <- response[,"time"]
    # unique event times
    times <- unique(times)
  }

  # }}}
  # {{{ newdata
  if(missing(newdata)){
    newdata <- eval(x$call$data)
  }
  ## stop("newdata argument is missing")

  # }}}
  # {{{ xlim, ylim

  if (missing(xlim)) xlim <- c(0, max(times))
  at <- times <= xlim[2]
  orig.X <- times[at]
  X <- times[at]
  
  # }}}  
  # {{{ predict newdata at times
  y <- predictEventProb(object=x,
                        newdata=newdata,
                        times=orig.X,
                        cause=cause)
  
  # }}}
  # {{{  plot arguments

  nlines <- NROW(y)
  
  if (missing(ylab)) ylab <- "Event probability"
  if (missing(xlab)) xlab <- "Time"
  if (missing(ylim)) ylim <- c(0, 1)
  if (missing(lwd)) lwd <- rep(3,nlines)
  if (missing(col)) col <- rep(1,nlines)
  if (missing(density)){
    if (nlines>5){
      density <- pmax(33,100-nlines)
    }
    else density <- 100
  }
  ## print(density)
  if (density<100){
    col <- sapply(col,function(coli){
      ccrgb=as.list(col2rgb(coli,alpha=TRUE))
      names(ccrgb) <- c("red","green","blue","alpha")
      ccrgb$alpha=density
      cc=do.call("rgb",c(ccrgb,list(max=255)))
    })
  }
  if (missing(lty)) lty <- rep(1, nlines)
  if (length(lwd) < nlines) lwd <- rep(lwd, nlines)
  if (length(lty) < nlines) lty <- rep(lty, nlines)
  if (length(col) < nlines) col <- rep(col, nlines)
  
  axis1.DefaultArgs <- list()
  
  axis2.DefaultArgs <- list(at=seq(0,1,.25))
         
  plot.DefaultArgs <- list(x=0,
                           y=0,
                           type = "n",
                           ylim = ylim,
                           xlim = xlim,
                           xlab = xlab,
                           ylab = ylab)
  
  
  
  legend.DefaultArgs <- list(legend=rownames(y),
                             lwd=lwd,
                             col=col,
                             lty=lty,
                             cex=1.5,
                             bty="n",
                             y.intersp=1.3,
                             x="topright")
  
  # }}}
  # {{{ smart controls

  if (match("legend.args",names(args),nomatch=FALSE)){
    legend.DefaultArgs <- c(args[[match("legend.args",names(args),nomatch=FALSE)]],legend.DefaultArgs)
    legend.DefaultArgs <- legend.DefaultArgs[!duplicated(names(legend.DefaultArgs))]
  }
  smartA <- prodlim:::SmartControl(call=list(...),
                                   keys=c("plot","legend","axis1","axis2"),
                                   ignore=c("x", "newdata", "times", "xlim","ylim","xlab","ylab","col","lty","lwd","add","legend","percent","axes","legend.args"),
                                   defaults=list("plot"=plot.DefaultArgs,
                                     "legend"= legend.DefaultArgs,
                                     "axis1"=axis1.DefaultArgs,
                                     "axis2"=axis2.DefaultArgs),
                                   forced=list("plot"=list(axes=FALSE),
                                     "axis1"=list(side=1),
                                     "axis2"=list(side=2)),
                                   verbose=TRUE)

  # }}} 
  # {{{ empty plot
  
  if (!add) {
    do.call("plot",smartA$plot)
    if (axes){
      do.call("axis",smartA$axis1)
      if (percent & is.null(smartA$axis1$labels))
        smartA$axis2$labels <- paste(100*smartA$axis2$at,"%")
      do.call("axis",smartA$axis2)
    }
  }

  # }}}
  # {{{ adding lines
  nix <- lapply(1:nlines, function(s) {
    lines(x = X, y = y[s,], type = "s", col = col[s], lty = lty[s], lwd = lwd[s])
  })

  # }}}
  # {{{ legend

  if(legend==TRUE && !add && !is.null(rownames(y))){
    save.xpd <- par()$xpd
    do.call("legend",smartA$legend)
    par(xpd=save.xpd)
  }
  
  # }}}

  invisible(x)
}

