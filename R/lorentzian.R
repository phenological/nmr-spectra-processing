#' Lorentzian function

#' Interpolates a lorentzian with the given parameters on the given points. 
#' Vectorized
#' #' @export
#' @param x numeric, points to be interpolated
#' @param Max numeric, maxima of the lorentzian(s)
#' @param Mean numeric, mean(s) of the lorentzian(s)
#' @param fwhm numeric, full-width-at-half-max(s) of the lorentzians(s)
#' @returns numeric matrix with lorentzian(x) for each combination of parameters
#' @examples lorentzian(1:20,Max=1,Mean=c(10,15))
lorentzian <- function(x,max=1,mean=0,fwhm=1){
  gamma2 <- fwhm/2
  gamma2 <- gamma2^2
  mapply(function(max.,mean.,gamma2.) max. * gamma2./((x-mean.)^2+gamma2.)
         , max, mean, gamma2)
}