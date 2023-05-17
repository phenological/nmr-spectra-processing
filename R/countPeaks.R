#' Counts the number of peaks in a data.frame of peaks
#' 
#' @export
#' qol wrapper to dim that returns 0 for failed peak picking (no
#' data.frame produced)
#' @param peaks data.frame of peaks. Run about.NMRPeak1D() for information on
#' the format of peaks
countPeaks <- function(peaks){
  if (is.data.frame(peaks)) dim(peaks)[1]
  else 0
}