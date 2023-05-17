#' Pseudo-Voigt function
#' 
#' Interpolates a Voigt function, using the pseudo-Voigt approximation, 
#' on the given points with the given parameters. 
#' The pseudo-Voigt approximation approximates a Voigt function as a linear 
#' combination of a lorentzian and a gaussian. 
#' Vectorized.
#' @export
#' @param x numeric points to be interpolated
#' @param mean numeric, mean of the distribution
#' @param fwhm numeric, full width at half max of the distribution
#' @param mu numeric, <= 1. The mu parameter of the pseudo-Voigt approximation.
#' The coefficientes of the linear combination are mu for the lorentzian and
#' 1 - mu for the gaussian
#' @returns numeric matrix with gaussian(x) for each combination of parameters
#' @examples pseudoVoigt(1:20,Max=1,Mean=c(10,15),mu=c(0,0.5))
pseudoVoigt <- function(x,max=1,mean=0,fwhm=1,mu=0){
  mapply(function(max.,mean.,fwhm.,mu.)
    mu.*lorentzian(x,max.,mean.,fwhm.) + (1-mu.)*gaussian(x,max.,mean.,fwhm.)
    , max, mean, fwhm, mu)
}