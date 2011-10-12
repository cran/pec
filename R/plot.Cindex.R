plot.Cindex <- function(x,ylim=c(.4,1),xlim=c(0,x$maxtime),abline=TRUE,...){
  argList <- match.call(expand=TRUE)
  argList[[1]] <- as.name("list")
  argList <- eval(argList,parent.frame())
  argList <- c(list("predErr"=switch(x$method$internal.name,"noPlan"={"AppCindex"},
                      paste(x$splitMethod$internal.name,"Cindex",sep=""),ylab="C-index")),argList)
  argList$ylim <- ylim
  argList$xlim <- xlim
  argList$x$exact <- FALSE
  do.call("plot.pec", argList)
  if (abline==TRUE)
    abline(h=.5,col="gray",lty=3,lwd=3,xpd=FALSE)
}
