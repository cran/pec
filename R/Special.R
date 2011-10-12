Special <-  function(x,
                     y,
                     addprederr,
                     models,
                     bench,
                     benchcol,
                     times,
                     maxboot,
                     bootcol,
                     col,
                     lty,
                     lwd){
   if(length(models) != 1)
     stop("Need to choose one and only one 'models'")
   at <- sindex(x$time,times)

   # add the bootstrap curves 
   boot <- x$BootstrapCrossValErrMat[[models]]
  if(is.null(maxboot)) maxboot <- NCOL(boot) 
   
  nix <- lapply(1:maxboot,function(b){
    lines(times,boot[b,at,drop=TRUE],col=bootcol,type="s",lwd=2.5)
  })

   # add the predErr of models
   addw <- lines(times, y[[models]][at], col=col[1], lty=lty[1], lwd=lwd[1])
   
   # add the xtra chosen whats of models
   if(!is.null(addprederr)){
     nx <- length(addprederr)
     addx <- lapply(1:nx, function(a){
       newy <- x[[addprederr[a]]][[models]]
       lines(times, newy[at], col=col[a+1], lwd=lwd[a+1], lty=lty[a+1])
     })
   }
   
   # add the benchmark models
   if(bench!= FALSE){
     addb <- lines(times, y=y[[bench]][at], col=benchcol, lwd=lwd[1])
   }
 }

