#' Aligns series to a reference
#' 
#' Computes the shift that maximizes each series' cross-correlation to the reference
#' @param x numeric vector or matrix, series in rows
#' @param ref logical, character, numeric, function. Specifies the reference for alignment. May be a series of the same length as the input series (only option if x is a vector), the index of the row of x to be used as reference, a logical filter that selects a unique row of x to be used as reference, a function to compute the reference from x, or a character representing an implemented reference. Valid character references are 'median' (default: reference is the median of x) and 'mean'.
#' @param threshold numeric, cross-correlation (similarity) must be above this threshold for the spectrum to be aligned, otherwise it is left unshifted.
#' @param shift logical. If TRUE, (default) returns the aligned series. If FALSE, returns the shifts that align the series
#' @param padding character, the method to be used to fill the empty extremes of the shifted series. See 'pad' for details. Default: "zeroes".
#' @param from numeric, optional. Filter for the region to be used for padding after the series is shifted in the "sampling" method. See 'pad' for details. Default: last 1/15th points
#' @param plot logical. Argument passed to ccf; if TRUE, each spectrum's correlation to the reference is plotted. Default: FALSE. You should not need to change the default unless you are getting misalignments and you want to check the cross-correlation to fine-tune the alignment parameters.
#' @param ... additional arguments to ccf, see Details
#' @returns either a numeric vector of shifts or a matrix with the shifted series in its rows
#' @export
alignSeries <- function(x, ref=c("median","mean","other options in the documentation")[1]
                        ,threshold=0.6, shift=TRUE
                        ,padding=c("zeroes","circular","sampling")[1]
                        ,from = as.integer(length(x)*14/15):length(x)
                        ,plot = FALSE, ...){
  alignSeries.numeric <- function(x, ref, threshold, shift,padding,from,plot
                                  , ...){
    cc <- ccf(x, ref, type="correlation", plot=plot, ...)
    ccmax <- which.max(cc$acf)
    if (length(ccmax)==0){
      cat(crayon::yellow("nmr.spectra.processing::alignSeries>>"
                         ,"Cross-correlation diverged,"
                         ,"may be due to non-converging input"
                         ,"(i.e. convex or constant series)."
                         ,"Series left unshifted.\n"))
      delta = 0
    }
    else{
      delta <- as.vector(cc$lag)[ccmax]
      if (cc$acf[ccmax] < threshold){
        cat(crayon::yellow("nmr.spectra.processing::alginSeries>>"
                           ,"Cross-correlation lower than threshold"
                           ,"Series left unshifted.\n"))
        delta = 0
      }
    }
    if (shift) return(shiftSeries(x,-delta,padding=padding,from=from))
    return(-delta)
  }
  
  if (is.vector(x)){
    if (!is.numeric(ref) | length(ref) != length(x)){
      cat(crayon::red("nmr.spectra.processing::alingSeries>>"
                      ,"expected a numeric reference of the same length"
                      ,"as the input series\n"))
      stop()
    }
    return(alignSeries.numeric(x,ref,threshold,shift,padding,from,plot,...))
  }
  
  #Parse ref and build reference spectrum
  # if(!(is.logical(ref)|is.character(ref)|is.function(ref)|is.numeric(ref))){
  #   cat(crayon::red("nmrSpectraProcessing::alignSeries >>", "invalid ref type\n"))
  #   stop()
  # }
  if (is.function(ref)){
    ref <- ref(x)
  }
  if (is.character(ref)){
    if (ref=="median") ref <- apply(x,2,median)
    else{
      if (ref=="mean") ref <- apply(x,2,mean)
      else {
        cat(crayon::red("nmrSpectraProcessing::alignSeries >>"
                        ,"Invalid character ref: must be 'median' (default)"
                        ,"or 'mean'\n"))
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
  return(apply(x,1,function(v){
    alignSeries.numeric(v, ref=ref, threshold=threshold, shift=shift
                        ,padding=padding,from=from, plot = plot, ...)
    })
  )
}
