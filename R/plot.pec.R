"plot.pec" <-  function(x,
                        what="PredErr",
                        who,
                        crps=FALSE,
                        xlim=c(x$start,x$maxtime),
                        ylim=c(0,0.3),
                        xlab="Time",
                        ylab,
                        lwd.lines=2,
                        axes=TRUE,
                        col,
                        lty,
                        lines.type,
                        smooth=FALSE,
                        add.refline=FALSE,
                        add=FALSE,
                        legend=ifelse(add,FALSE,TRUE),
                        legend.text,
                        legend.args=list(),
                        specials=NULL,...){
  if (missing(lines.type)) if (!x$exact || smooth) lines.type <- "l" else lines.type <- "s"
  at <- x$time <= xlim[2]
  X <- x$time[at]
  if (missing(who)) who <- 1:length(x$models)
  else
    if (!is.numeric(who))
      who <- names(x$models)[match(who,names(x$models))]
  if (what=="special"){
    if(!(x$method$internal.name=="boot632plus"||x$method$internal.name=="boot632"))
      stop("Plotting method `special' requires prediction error method `boot632plus' or `boot632'")
    if (is.null(specials$bootcol)) specials$bootcol <- "gray77"
    mc <- match.call()
    mc$legend <- FALSE
    if (is.null(specials$what)) specials$what <- c("PredErr")
    if (!is.null(specials$bench)){
      if (is.null(specials$benchcol)) specials$benchcol <- 1
      who <- who[who!=specials$bench]
      mcbench <- mc
      mcbench$who <- specials$bench
      mcbench$what <- "PredErr"
      mcbench$add <- FALSE
    }
    lapply(who,function(w){
      mcw <- mc
      mcw$who <- w
      if (!is.null(specials$bench)){
        eval(mcbench)
      }
      else{
        mcw$add <- FALSE
        mcw$what <- "PredErr"
        mcw$lty <- if (is.null(specials$lty.632Err)) 1 else specials$lty.632Err
        eval(mcw)
      }
      boot <- x$OutOfBagErrMat[[w]]
      nix <- lapply(1:NCOL(boot),function(b){
        lines(X,boot[,b,drop=TRUE],col=specials$bootcol,type="s",lwd=2.5)
      })
      mcw$col <- if (is.null(specials$col)) 2 else specials$col
      mcw$add <- TRUE
      if ("AppErr" %in% specials$what){
        mcw$what <- "AppErr"
        mcw$lty <- if (is.null(specials$lty.AppErr)) 2 else specials$lty.AppErr
        eval(mcw)
      }
      if ("OutOfBagErr" %in% specials$what){
        mcw$what <- "OutOfBagErr"
        mcw$lty <- if (is.null(specials$lty.OutOfBagErr)) 3 else specials$lty.OutOfBagErr
        eval(mcw)
      }
      if (!is.null(specials$bench)){
        mcbench$add <- TRUE
        eval(mcbench)
      }
      mcw$what <- "PredErr"
      mcw$lty <- if (is.null(specials$lty.632Err)) 1 else specials$lty.632Err
      eval(mcw)
    })
    invisible(x)
  }
  else{
    if (crps==TRUE)
      y <- t(crps(x,who=who,what=what,times=X,print=FALSE))
    else
      y <- do.call("cbind",x[[what]][who])[at,,drop=FALSE]
    ## y <- x[[what]][at,who,drop=FALSE]

    if (length(y)==0) stop("No plotting values: check if x[[what]][who] is a list of numeric vectors.")
    
    if (!add)
      if (missing(ylab)) ylab <- switch(grep(what,c("PredErr","R2","overfit","AppErr","OutOfBagErr")),
                                        "Prediction error",
                                        expression(R^2),
                                        "Relative overfit",
                                        "Prediction error",
                                        "Prediction error")
    uyps <- unlist(y)
    uyps <- uyps[!is.infinite(uyps)]
    max.y <- max(uyps,na.rm=T)
    ymax <- max(max.y,ylim[2])
    xmin <- min(X)
    if (max.y>ylim[2])
      ylim <- if (what=="PredErr")
        c(0,ceiling(ymax*10)/10)
      else
        c(0,ceiling(max(unlist(y),na.rm=T)*10))/10
    ##   ceiling(c(min(unlist(y),na.rm=T),max(unlist(y),na.rm=T)*10))/10
    if (!add){
      plot(0,0,ylab=ylab,xlab=xlab,type="n",ylim=ylim,xlim=xlim,axes=FALSE,...)
      if (axes==TRUE){
        ##       if (what=="pred.error" || max(unlist(y)>.25,na.rm=T))
        ##         axis(2,at=seq(0,0.25,.05),labels=seq(0,0.25,.05))
        ##       else
        axis(2)  
        axis(1)
      }
    }
    
    # adding the curves 
    # --------------------------------------------------------------------

    nfit <- ncol(y)
    if (missing(lty)) lty <- 1
    if(length(lty)<nfit) lty <- cbind(lty,1:nfit)[,1]
    if (missing(col)) col <- 1:nfit
    if(length(col)<nfit) col <- cbind(col,1:nfit)[,1]
    for (f in 1:ncol(y)){
      Yvals <- y[,f]
      if (smooth==TRUE){
        Yvals[!is.na(Yvals) & !is.infinite(Yvals)] <- smooth(Yvals[!is.na(Yvals) & !is.infinite(Yvals)],kind="3R")
      }
      points(X,Yvals,type=lines.type,col=col[f],lty=lty[f],lwd=lwd.lines)
    }

    if (add.refline) abline(h=.25,lty=3,lwd=2,col=1)
  
    # legend
    # --------------------------------------------------------------------

    if (legend){
      old.xpd <- par()$xpd
      par(xpd=TRUE)
      if (missing(legend.text) && !("legend" %in% names(legend.args)))
        legend.args <- c(list(legend=names(x$models)[who]),legend.args)
      else{
        if (!("legend" %in% names(legend.args)) && length(legend.text==length(x$models)))
          legend.args <- c(list(legend.text),legend.args)
      }
      
      if (!("x" %in% names(legend.args))) legend.args <- c(list(x=0),legend.args)
      if (!("y" %in% names(legend.args))) legend.args <- c(list(y=(ylim+.1*ylim)[2]),legend.args)
      if (!("col" %in% names(legend.args))) legend.args <- c(list(col=col),legend.args)
      if (!("lty" %in% names(legend.args))) legend.args <- c(list(lty=lty),legend.args)
      if (!("lwd" %in% names(legend.args))) legend.args <- c(list(lwd=lwd.lines),legend.args)
      if (!("cex" %in% names(legend.args))) legend.args <- c(list(cex=1.5),legend.args)
      if (!("bty" %in% names(legend.args))) legend.args <- c(list(bty="n"),legend.args)

      do.call("legend",legend.args)  
      invisible(x)
      par(xpd=old.xpd)
    }
  }
}
