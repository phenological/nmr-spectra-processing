#' Baseline Correction
#' 
#' Wrapper to \code{\link[ptw]{baseline.corr}}. Estimates baseline by asymmetric
#'  least squares and subtracts it from the spectra
#' @param y numeric or matrix, single spectrum intensities or intensities matrix
#' with spectra in rows
#' @param ... other arguments to \code{\link[ptw]{asysm}}
#' @details Presently this is \code{\link[ptw]{baseline.corr}} straight out of
#'  the box, wrapped and re-packaged for convenience.
#' @returns a vector or matrix of the same dimensions as the input with 
#' baseline-corrected spectra.
#' @importFrom ptw baseline.corr
#' @export
baselineCorrection <- function(y, ...) {
  bc <- baseline.corr(y,...)
  return(bc)
}
