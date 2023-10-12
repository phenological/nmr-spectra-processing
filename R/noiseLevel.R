#' Estimates the noise level of the spectrum
#' 
#' @param ppm numeric, chemical shift scale of the spectrum.
#' @param y numeric, spectrum intensities.
#' @param level numeric, quantile to estimate noise, specified as a probability
#' as in \code{\link[stats]{quantile}}. Default: 0.99.
#' @param rOref numeric, limits of the Region of Reference used to estimate 
#' noise. Default: 9.8 - 10 ppm.
#' @details Takes a 'blank' region of the spectrum as reference and
#'  estimates noise as the \code{level}-probability quantile of intensities
#'  in that region.
#' @returns numeric, estimation of the spectrum's level of noise.
#' @importFrom stats quantile
#' @export
noiseLevel <- function(ppm,y,level=0.99,rOref=c(9.8,10)){
  fi <- ppm >= rOref[1] & ppm <= rOref[2]
  if (is.matrix(y)){
    return(apply(y[,fi],1,function(r){
      quantile(r,level,names=FALSE)
    })
    )
  }
  return(quantile(y[fi],level,names=FALSE))
}