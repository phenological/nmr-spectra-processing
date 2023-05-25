#' Baseline Correction
#' @param Y matrix, series in rows
#' @param ... other arguments to ptw::asysm
#' @returns a matrix of same size with each row being corrected
#' @importFrom ptw baseline.corr
#' @details ptw:baseline.corr wrapper. Estimates baseline by asymmetric least squares using ptw::asysm and subtracts it from the spectra
#' @export
baselineCorrection <- function(Y, ...) {
  bc <- ptw::baseline.corr(Y,...)
  return(bc)
}
