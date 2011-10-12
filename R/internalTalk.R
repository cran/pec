internalTalk <- function(x,y,sign="'"){
  if (y>100){
    if (y<500){
      if (x %in% seq(0,y,10)) cat(paste("\n",x))
    }
    else{
      if (x %in% seq(0,y,100)) cat(paste("\n",x))
    }
  }
  else{
    if (y>10){
      if (x %in% seq(0,y,10)) cat(paste("\n",x))
    }
    else
      if (y>1)
        cat(paste(sign))
  }
}

