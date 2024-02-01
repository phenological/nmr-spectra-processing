#' Align series to a reference
#' 
#' Invokes \code{\link[stats]{ccf}} to compute the cross-correlations
#' between the series to be aligned and the reference series and extracts the
#'  shifts that maximize the cross-correlations
#' @param x numeric vector or matrix, series in rows.
#' @param ref logical, character, numeric or function. Specifies the reference
#'  series for alignment. It is interpreted according to type:\itemize{
#'  \item Numeric vector of the same length as the input series: literal reference
#'  \item A number \code{<= dim(x)[1]}: row index of the reference series
#'  \item logical of length \code{== dim(x)[2])}: logical filter that selects a
#'  single row of \code{x} to be used as reference
#'  \item function: the function is called with \code{x} as argument, the 
#'  returned value is used as reference
#'  \item character: "median" (default) or "mean", use the median ir mean series
#'   as reference.}
#' @param threshold numeric, the cross-correlation must be above this
#'  threshold for the series to be aligned, otherwise it is left unshifted.
#'  Increase this value if you know there are very big shifts. Decrease it if 
#'  you observe over-shifting.
#' @param shift logical. If TRUE (default), returns the aligned series. If FALSE,
#'  returns the shifts that align the series.
#' @param padding character, the method to be used to fill the empty extremes of
#'  the shifted series. See \code{\link{pad}} for details. Default: "zeroes".
#' @param from numeric, optional. Filter for the region to be used for padding
#'  after the series is shifted in the "sampling" method. See \code{\link{pad}}
#'  for details. Default: last 1/15th points.
#' @param ... additional arguments to \code{\link[stats]{ccf}}
#' @returns either a numeric vector of shifts or a matrix with the shifted series
#'  in its rows
#'  @import methods
#' @importFrom stats ccf median
#' @export
alignSeries <- function(x, ref=c("median","mean","more options in documentation")[1]
                        ,threshold=0.6, shift=TRUE
                        ,padding=c("zeroes","circular","sampling")[1]
                        ,from = as.integer(length(x)*14/15):length(x), ...){
  alignSeries.numeric <- function(x, ref, threshold, shift,padding,from, ...){
    # if(!("plot" %in% names(list(...))))
    cc <- ccf(x, ref, type="correlation", plot=FALSE, ...)
    ccmax <- which.max(cc$acf)
    if (length(ccmax)==0){
      cat(crayon::yellow("nmr.spectra.processing::alignSeries >>"
                         ,"Cross-correlation diverged,"
                         ,"may be due to non-converging input"
                         ,"(i.e. convex or constant series)."
                         ,"Series left unshifted.\n"))
      delta = 0
    }
    else{
      delta <- as.vector(cc$lag)[ccmax]
      if (cc$acf[ccmax] < threshold){
        cat(crayon::yellow("nmr.spectra.processing::alginSeries >>"
                           ,"Cross-correlation", cc$acf[ccmax]
                           ,"lower than threshold", threshold
                           ,".\nSeries left unshifted.\n"))
        delta = 0
      }
    }
    if (shift) return(shiftSeries(x,-delta,padding=padding,from=from))
    return(-delta)
  }
  
  if (is.vector(x)){
    if(!is.numeric(x)){
      cat(crayon::yellow("nmr.spectra.processing:alignSeries >>"
                         ,"Non-numeric argument x being cast as.numeric\n"
                         ,"Unpredictable results will follow if casting to"
                         ,"numeric vector fails\n"))
      x <- as.numeric(x)
    }
    if (!is.numeric(ref)){
      cat(crayon::yellow("nmr.spectra.processing:alignSeries >>"
                         ,"Exoected numeric reference\n"
                         ,"Non-numeric argument ref being cast as.numeric\n"
                         ,"Unpredictable results may follow if casting to"
                         ,"numeric vector fails\n"))
      ref <- as.numeric(ref)
    }
    if (length(ref) != length(x)){
      cat(crayon::red("nmr.spectra.processing::alingSeries >>"
                      ,"expected a numeric reference of the same length"
                      ,"as the input series\n"))
      stop()
    }
    return(alignSeries.numeric(x,ref,threshold,shift,padding,from,...))
  }
  
  #Parse input as matrix
  if (!is.matrix(x)){
    cat(crayon::yellow("nmr.spectra.processing:alignSeries >>"
                       ,"Argument x being cast as.matrix\n"
                       ,"Unpredictable results may follow if casting to"
                       ,"numeric matrix fails\n"))
    x <- as.matrix(x)
  }
  if (!is.numeric(x)){
    cat(crayon::red("nmr.spectra.processing::alingSeries >>"
                    ,"expected x to be a numeric vector or matrix \n"))
    stop()
  }
  
  #Parse ref and build reference spectrum
  if (is.function(ref)){
    ref <- ref(x)
  }
  if (is.character(ref)){
    if (ref=="median") ref <- apply(x,2,median)
    else{
      if (ref=="mean") ref <- apply(x,2,mean)
      else {
        cat(crayon::red("nmrSpectraProcessing::alignSeries >>"
                        ,"Invalid character ref: must be 'median' (default),"
                        ,"'mean', or a function.\n"))
        stop()
      }
    }
  }
  
  if (is.numeric(ref)){
    if (length(ref) == 1){
      if (ref <= dim(x)[1]) ref <- x[ref,]
      else{
        cat(crayon::red("nmrSpectraProcessing::alignSeries >>"
                        ,"ref index not found in input matrix\n"))
      }
    }
  }
  
  if (is.logical(ref)) ref <- x[ref,]
    
  if (!is.numeric(ref) | length(ref) != dim(x)[2]) {
    cat(crayon::red("nmrSpectraProcessing::alignSeries >>",
                        "Invalid reference\n"))
    stop()
  }
  return(t(
    apply(x,1,function(v){
      alignSeries.numeric(v, ref=ref, threshold=threshold, shift=shift
                          ,padding=padding,from=from, ...)
      })
    )
  )
}
