#' Fast spectra overlay plot
#' 
#' Wrapper to matplot for plotting spectra. Crude compared to ggplot2 but a lot faster.
#' @param x numeric, your x scale e.g. ppm
#' @param Y matrix, spectra matrix, one spectrum per row
#' @param limits, numeric, optional, limits of the region to be plotted. Defaults to full spectra.
#' @param by, numeric, optional, number of spectra to be overlay on each plot
#' @param type, optional, defaults to "l" (lines), default strongly recommended
#' @param ..., additional options to be passed to matplot
#' @returns NULL
#' @export
#Candidate to utilities.R
smatplot <- function(x, Y, by, limits, type="l",...){
  if (!missing(limits)){
    fi <- x >= limits[1] & x <= limits[2]
    x <- x[fi]
    Y <- Y[,fi]
  }
  if (missing(by)) matplot(x,t(Y),type=type,...)
  else{
    n <- floor(dim(Y)[1] / by)
    r <- dim(Y)[1] %% by
    for (j in 1:n) matplot(x,t(Y[1:by + by*(j-1),]),type=type,...)
    if (r>0) matplot(x,t(Y[(by*n+1):(by*n+r),]),type=type,...)
  }
}