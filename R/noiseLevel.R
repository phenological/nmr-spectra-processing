#' Estimates noise level of the spectrum
#' 
#' @param x numeric, your x scale e.g. ppm
#' @param y numeric, spectrum intensities
#' @param level numeric, quantile to estimate noise, see details
#' @param rOref numeric, limits of the region Of reference to be used to estimate noise, see details
#' @details Takes a 'blank' region of the spectrum as reference (rOref) and estimates noise as the level% quantile of intensity points in that region
#' @returns numeric, estimation of the spectrum's level of noise, namely the level% quantile of intensities in the rOref
#' @export
noiseLevel <- function(x,y,level=0.99,rOref=c(9.8,10)){
  fi <- x >= rOref[1] & x <= rOref[2]
  quantile(y[fi],level,names=FALSE)
}