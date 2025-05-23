#' Apply a phase correction to a complex NMR spectrum. phi1 is specified at index=0
#' 
#' @param yr list, spectrum intensities, real part
#' @param yi list, spectrum intensities, imaginary part
#' @param phi0 numeric, phase 0 order correction.
#' @param phi1 numeric, phase 1 order correction.
#' @param reverse boolean, TRUE if the chemical shift scale is inverted 
#' (decreases from left to right). Default=FALSE
#' @returns list, A list with 2 elements corresponding to the corrected spectrum (yr, yi)
#' @export
phaseCorrection <- function(yr, yi, phi0, phi1, reverse=FALSE) {
  phi0 <- phi0 * pi / 180
  phi1 <- phi1 * pi / 180
  
  lng <- length(yr)
  delta <- phi1 / lng
  firstAngle <- phi0
  
  if (!reverse) {
    delta <- -delta;
    firstAngle <- firstAngle + phi1;
  }
  
  
  alpha <- 2 * sin(delta / 2) ** 2;
  beta <- sin(delta)
  cosTheta <- cos(firstAngle);
  sinTheta <- sin(firstAngle);
  
  yr2 <- list();
  yi2 <- list();
  for (i in 1:lng) {
    yr2[i] <- yr[i] * cosTheta - yi[i] * sinTheta;
    yi2[i] <- yi[i] * cosTheta + yr[i] * sinTheta;
    
    tmpCos <- cosTheta - (alpha * cosTheta + beta * sinTheta)
    tmpSin <- sinTheta - (alpha * sinTheta - beta * cosTheta)
    cosTheta <- tmpCos;
    sinTheta <- tmpSin;
  }
  
  
  return(list(unlist(yr2), unlist(yi2)))
}