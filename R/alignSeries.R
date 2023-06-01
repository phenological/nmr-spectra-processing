#' Aligns series to a reference
#' 
#' Computes the shift that maximizes each series' cross-correlation to the
#' @param Y matrix, series in rows
#' @param ref character, numeric or function. Specifies the reference for alignment. May be a series of the same length as the input series, the index of the row of Y to be used as reference, a function to compute the reference from the input, or a character representing an implemented reference. Valid character references are 'median' (default: reference is the median of the input) and 'mean'. If ref is a function, it run on Y and the result is used as reference.
#' @param threshold numeric, cross-correlation (similarity) must be above this threshold for the spectrum to be aligned; otherwise, it is left unshifted.
#' @param shift logical. If TRUE, (default) returns the aligned series. If FALSE, returns the shifts that align the series
#' @param method character, the method to be used to fill the empty extremes of the shifted series (padding). See 'pad' for details. Default: "zeroes".
#' @param using numeric, the number of points from the extremes to be used for padding in the "sampling" method. See 'pad' for details. Default: 1/15th of the series' length.
#' @param plot logical. Argument passed to ccf; if TRUE, each spectrum's correlation to the reference is plotted. Default: FALSE. You should not need to change the default unless you are getting misalignments and you want to check the cross-correlation to fine-tune the alignment parameters.
#' @param ... additional arguments to ccf, see Details
#' @returns either a numeric vector of shifts or a matrix with the shifted series in its rows
#' @export
alignSeries <- function(Y, ref="median", threshold=0.6, shift=TRUE
                        , method="zeroes", using=(dim(Y)[2]*14/15:dim(Y)[2])
                        , from=apply(Y,2,median), plot = FALSE, ...){
  #QC input Y
  if(!is.matrix(Y)){
    cat(crayon::red("alignSeries >>","Y must be a matrix\n"))
    stop()
  }
  #Parse ref and build reference spectrum
  if(!(is.logical(ref)|is.character(ref)|is.function(ref)|is.numeric(ref))){
    cat(crayon::red("alignSeries >>", "invalid ref type\n"))
    stop()
  }
  if (is.function(ref)){
    ref <- ref(Y)
  }
  if (is.character(ref)){
    if (ref=="median") ref <- apply(Y,2,median)
    else{
      if (ref=="mean") ref <- apply(Y,2,mean)
      else {
        cat(crayon::red("alignSeries >>",
                        "Invalid character ref: must be 'median; (default) or 'mean'\n"))
        stop()
      }
    }
  }
  if (is.numeric(ref)){
    if (length(ref) == 1){
      if (ref <= dim(Y)[1]) ref <- Y[ref,]
      else{
        cat(crayon::red("alignSeries >>"
                        ,"ref index not found in input matrix\n"))
      }
    }
  }
    
  if (length(ref) != dim(Y)[2]) {
    cat(crayon::red("alignSeries >>",
                        "Invalid ref: reference spectrum does not match the length of input spectra\n"))
    stop()
  }
  #print(length(ref))
  #print(ref[1:10])
  t(apply(Y,1,function(y){
    cc <- ccf(y, ref, type="correlation", plot=plot, ...)
    ccmax <- which.max(cc$acf)
    delta <- as.vector(cc$lag)[ccmax]
    if (cc$acf[ccmax] < threshold) delta = 0
    if (shift) return(shiftSeries(y,-delta,method=method,using=using,from=from))
    return(-delta)
  }))
}
