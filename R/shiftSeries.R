#' 
shiftSeries <- function(x,shift,padding="sampling"
                        ,from=x,using=as.integer(length(x)*14/15):length(from)
){
  direction <- sign(shift)
  shift <- abs(shift)
  padded <- pad(x,shift,-direction,method=padding,using=using,from=from)
  if (direction==1) return(padded[1:length(x)])
  return(padded[(shift+1):length(padded)])
}