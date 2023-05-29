shiftSeries <- function(x,shift,method="sampling",using=length(x)/15){
  direction <- sign(shift)
  shift <- abs(shift)
  padded <- pad(x,shift,-direction,method,using)
  if (direction==1) return(padded[1:length(x)])
  return(padded[(shift+1):length(padded)])
  #trail_l <- if (shift <= 0) 0 else shift
  #trail_r <- if (shift < 0) -shift else 0
  #c(rep(0, trail_l), y[(trail_r+1):(length(y)-trail_l)], rep(0,trail_r))
}