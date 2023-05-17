#' Baseline Correction
#' @param Y matrix, series in rows
#' @param ... other arguments to ptw::asysm
#' @returns a matrix of same size with each row being corrected
#' @importFrom ptw baseline.corr
#' @export
#' @details ptw:baseline.corr wrapper. Estimates baseline by asymmetric least squares using ptw::asysm an subtracts it from the spectra
baselineCorrection <- function(Y, ...) {
  bc <- ptw::baseline.corr(Y,...)
  return(bc)
}
