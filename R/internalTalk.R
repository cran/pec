internalTalk <- function(x,y){
  if (y>100){
    if (x %in% seq(0,y,100)) cat(paste("\n",x))
  }
  else{
    if (y>10){
      if (x %in% seq(0,y,10)) cat(paste("\n",x))
    }
    else
      if (y>1)
        cat(paste("'"))
  }
}

