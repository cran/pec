# CFOREST
# --------------------------------------------------------------------


#' S3-wrapper function for cforest from the party package
#' 
#' S3-wrapper function for cforest from the party package
#' 
#' See \code{cforest} of the \code{party} package.
#' 
#' @param formula See \code{cforest} of the \code{party} package
#' @param data See \code{cforest} of the \code{party} package
#' @param \dots See \code{cforest} of the \code{party} package
#' @references Ulla B. Mogensen, Hemant Ishwaran, Thomas A. Gerds (2012).
#' Evaluating Random Forests for Survival Analysis Using Prediction Error
#' Curves. Journal of Statistical Software, 50(11), 1-23. URL
#' http://www.jstatsoft.org/v50/i11/.
#' @keywords survival
#' @export pecCforest
pecCforest <- function(formula,data,...){
    require(party)
    out <- list(forest=party::cforest(formula,data,...))
    class(out) <- "pecCforest"
    out$call <- match.call()
    out  
}


##' @S3method predictSurvProb pecCforest
predictSurvProb.pecCforest <- function (object, newdata, times, ...) {
    require(party)
    survObj <- party::treeresponse(object$forest,newdata=newdata)
    p <- do.call("rbind",lapply(survObj,function(x){
        predictSurvProb(x,newdata=newdata[1,,drop=FALSE],times=times)
    }))
    if (NROW(p) != NROW(newdata) || NCOL(p) != length(times)) 
        stop(paste("\nPrediction matrix has wrong dimension:\nRequested newdata x times: ",NROW(newdata)," x ",length(times),"\nProvided prediction matrix: ",NROW(p)," x ",NCOL(p),"\n\n",sep=""))
    p
}
