penalizedS3 <- function(formula,data,...){
  # {{{ distangle the formula
  ff <- as.character(formula)
  response <- formula(paste(ff[[2]],"~1",sep=""))
  terms <- strsplit(ff[[3]],"\\+|\\-")[[1]]
  terms <- sapply(terms,function(tt){## remove whitespace
    gsub(" ","",tt)
  })
  strippedTerms <- strsplit(terms,"[()]")
  # }}}
  # {{{ extract the penalized and unpenalized parts
  penalTerms <- sapply(strippedTerms,function(x){length(x)==2 && x[[1]]=="pen"})
  unpenalVarnames <- unlist(strippedTerms[penalTerms==FALSE])
  if (length(unpenalVarnames)>0){
    unpenalized <- formula(paste("~",paste(unpenalVarnames,collapse="+")))
    response <- update.formula(response,unpenalized)
  }
  penalizedVarnames <- unlist(sapply(strippedTerms[penalTerms==TRUE],
                                     function(x){strsplit(x[[2]],",")}),use.names=FALSE)
  penalizedVarPositions <- unlist(lapply(penalizedVarnames,function(x){
    if (length(splitter <- strsplit(x,":")[[1]])>1)
      seq(as.numeric(splitter)[1],as.numeric(splitter)[2],1)
    else
      match(x,names(data),nomatch=0)
  }),use.names=FALSE)
  penalizedVarPositions <- unique(penalizedVarPositions)
  ## print(penalizedVarPositions)
  if (any(tested <- (penalizedVarPositions>NCOL(data))|penalizedVarPositions<0))
    stop("Cannot find variable(s): ",names(data[tested]))
  penalized <- data[,penalizedVarPositions]
  # }}}
  # {{{ call S4 method
  fitS4 <- penalized(response=response,
                     penalized=penalized,
                     ## unpenalized=unpenalized,
                     data=data,
                     ...)
  # }}}
  # {{{ deliver S3 object
  fit <- list(fitS4=fitS4,call=match.call())
  class(fit) <- "penfitS3"
  fit
  # }}}
}
predictSurvProb.penfitS3 <- function(object,
                                     newdata,
                                     times,
                                     train.data,
                                     ...){
  penfit <- object$fit
  pCovaNames <- names(penfit@penalized)
  newPen <- newdata[,pCovaNames]
  ptemp <- predict(penfit,penalized=newPen,data=newdata)
  require(prodlim)
  pos <- sindex(jump.times=ptemp@time,eval.times=times)
  ## Remark: currently it is possible, but theoretically
  ## not allowed to carry predictions forward beyond the
  ## last jump.time
  p <- cbind(1,ptemp@curves)[,c(pos+1)]
  p
}

