# CTREE
pecCtree <- function(...){
 out <- list(ctree=party::ctree(...))
 class(out) <- "pecCtree"
 out$call <- match.call()
 out  
}

##' @S3method predictSurvProb pecCtree
predictSurvProb.pecCtree <- function (object, newdata, times, ...) {
    require(party)
    N <- NROW(newdata)
    NT <- length(times)
    survObj <- party::treeresponse(object$ctree, newdata=newdata)
    p <- do.call("rbind", lapply(survObj,function(x){
        predictSurvProb(x, newdata=newdata[1,,drop=FALSE], times=times)
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times))
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}
