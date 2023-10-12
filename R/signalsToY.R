#' Interpolates a \code{\linkS4class{NMRSignal1D}} on the given set of points
#' 
#' @param signal, \code{\linkS4class{NMRSignal1D}}, signal to interpolate
#' @param ppm, chemical shift values to interpolate
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
# importFrom nmr.peak.fitting gaussian lorentzian pseudoVoigt
#' @returns numeric vector of interpolated intensities
#' @export
signalToY <- function(signal,ppm){
  if("NMRSignal1D" %in% is(signal)){
    y <- 0
    fwhm <- signal@shape$params$fwhm
    mu <- signal@shape$params$mu
    for (p in signal@peaks)
      y <- y + pseudoVoigt(ppm, fwhm=fwhm, mu=mu, mean=p@x, max=p@y)
    return(y)  
  }
  cat(crayon::red("nmrSpectraProcessing::signalToY >>",
                  "Argument is not a NMRSignal1D object\n"))
}