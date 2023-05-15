#' Baseline Correction
#' @param Y matrix, series in rows
#' @returns a matrix of same size with each row being corrected
#' @importFrom ptw baseline.corr
baselineCorrection <- function(Y) {
  bc <- ptw::baseline.corr(Y)
  return(bc)
}
