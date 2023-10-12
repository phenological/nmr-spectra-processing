scaleFragment <- function(x,y,scale=c(2,2)) {
  if (length(x) != length(y)){
   cat(crayon::red("nmr-spectra-processing::scaleFragment >>"
                   ,"lengths of x and y do not match")) 
    stop()
  }
  #scale y
  ys <- y*scale[2]
  
  #scale x symmetrically with respect to the mid
  rangex <- range(x)
  dRange <- diff(range(x)) * (scale[1] - 1)
  newRange <- c(rangex[1] - dRange/2,rangex[2] + dRange/2)
  xs <- seq(newRange[1],newRange[2],length.out=length(x))
  
  #new scale for scaled x that is congruent with the original scale
  midx <- x[which.min(abs(x-mean(x)))]
  dx <- diff(x)[1]
  newx <- c(sort(seq(midx,newRange[1],by=-dx)),seq(midx+dx,newRange[2],by=dx))

  list(x=newx,y=interp1(xs,ys,newx,method="spline"))
}