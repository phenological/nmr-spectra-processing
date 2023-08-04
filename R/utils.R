#Initialization of the global v8 context within the package namespace
nsp.env <- new.env()
nsp.env$ct <- V8::v8()

#' Create a filter to crop a spectral region
#' 
#' Quality of life function to spare you a few uncomfortable key strokes. Questionable value.
#' @param x numeric, x scale e.g. ppm
#' @param roi numeric, optional. Upper and lower limit of the Region Of Interest (roi) to be cropped.
#' @param start, numeric, optional. Lower limit of the roi
#' @param end, numeric, optional. Upper limit of the roi
#' @details Either start, end or roi is required. Argument "roi" has priority. If start is given but not end,
#'  the lower limit is effectively set to min(x). If end is given but no start, the upper limit is effectively
#'  set to max(x).
#' @return logic, TRUE for elements of x within roi, FALSE otherwise

#Candidate to utilities.R
crop <- function(x,start=-Inf,end=Inf,roi){
  if (missing(roi))
    return(x >= start & x <= end)
  return(x >= roi[1] & x <= roi[2])
}