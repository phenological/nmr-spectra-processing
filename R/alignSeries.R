#' Aligns series to a reference
#' 
#' Computes the shift that maximizes each series' cross-correlation to the reference
#' @param y numeric vector or matrix, series in rows
#' @param ref character, numeric or function. Specifies the reference for alignment. May be a series of the same length as the input series, the index of the row of y to be used as reference, a function to compute the reference from the input, or a character representing an implemented reference. Valid character references are 'median' (default: reference is the median of the input) and 'mean'. If ref is a function, it run on y and the result is used as reference.
#' @param threshold numeric, cross-correlation (similarity) must be above this threshold for the spectrum to be aligned; otherwise, it is left unshifted.
#' @param shift logical. If TRUE, (default) returns the aligned series. If FALSE, returns the shifts that align the series
#' @param method character, the method to be used to fill the empty extremes of the shifted series (padding). See 'pad' for details. Default: "zeroes".
#' @param using numeric, the number of points from the extremes to be used for padding in the "sampling" method. See 'pad' for details. Default: 1/15th of the series' length.
#' @param plot logical. Argument passed to ccf; if TRUE, each spectrum's correlation to the reference is plotted. Default: FALSE. You should not need to change the default unless you are getting misalignments and you want to check the cross-correlation to fine-tune the alignment parameters.
#' @param ... additional arguments to ccf, see Details
#' @returns either a numeric vector of shifts or a matrix with the shifted series in its rows
#' @export
# alignSeries <- function(y,...){
#   UseMethod("alignSeries",y)
# }

alignSeries <- function(y, ref="median", threshold=0.6, shift=TRUE
                        ,padding="zeroes", using#=(dim(y)[2]*14/15:dim(y)[2])
                        ,from#=apply(y,2,median)
                        ,plot = FALSE, ...){
  alignSeries.numeric <- function(y, ref, threshold, shift
                                  ,padding="zeroes"
                                  ,from#=y
                                  ,using#=as.integer(length(y)*14/15):length(from)
                                  ,plot=FALSE, ...){
    
    cc <- ccf(y, ref, type="correlation", plot=plot, ...)
    ccmax <- which.max(cc$acf)
    if (length(ccmax)==0){
      cat(crayon::yellow("nmrSpectraProcessing::alignSeries>>","Cross-correlation diverged,"
                         ,"may be due to non-converging input (i.e. convex or constant series)."
                         ,"Series left unshifted."))
      delta = 0
    }
    else{
      delta <- as.vector(cc$lag)[ccmax]
      if (cc$acf[ccmax] < threshold) delta = 0
      if (shift) return(shiftSeries(y,-delta,padding=padding,using=using,from=from))
    }
    return(-delta)
  }
  
  if (is.vector(y)){
    if(missing(from))
      from <- y
    if(missing(using))
      using <- as.integer(length(y)*14/15):length(from)
    return(alignSeries.numeric(y,ref,threshold,shift,padding,from,using,plot,...))
  }
  
  #Parse ref and build reference spectrum
  if(!(is.logical(ref)|is.character(ref)|is.function(ref)|is.numeric(ref))){
    cat(crayon::red("nmrSpectraProcessing::alignSeries >>", "invalid ref type\n"))
    stop()
  }
  if (is.function(ref)){
    ref <- ref(y)
  }
  if (is.character(ref)){
    if (ref=="median") ref <- apply(y,2,median)
    else{
      if (ref=="mean") ref <- apply(y,2,mean)
      else {
        cat(crayon::red("nmrSpectraProcessing::alignSeries >>",
                        "Invalid character ref: must be 'median; (default) or 'mean'\n"))
        stop()
      }
    }
  }
  if (is.numeric(ref)){
    if (length(ref) == 1){
      if (ref <= dim(y)[1]) ref <- y[ref,]
      else{
        cat(crayon::red("nmrSpectraProcessing::alignSeries >>"
                        ,"ref index not found in input matrix\n"))
      }
    }
  }
    
  if (length(ref) != dim(y)[2]) {
    cat(crayon::red("nmrSpectraProcessing::alignSeries >>",
                        "Invalid ref: reference spectrum does not match the length of input spectra\n"))
    stop()
  }
  if (missing(from))
    from <- ref
  if (missing(using))
    using <- as.integer(dim(y)[2]*14/15):dim(y)[2]
  return(apply(y,1,function(v){
    alignSeries.numeric(v, ref=ref, threshold=threshold, shift=shift
                        ,padding=padding, using=using
                        ,from=from, plot = plot, ...)
    })
  )
  
  # if (is.matrix(res))
  #   return(t(res))
  # return(res)
  
    # cc <- ccf(y, ref, type="correlation", plot=plot, ...)
    # 
    # ccmax <- which.max(cc$acf)
    # if (length(ccmax)==0){
    #  cat(crayon::yellow("nmrSpectraProcessing::alignSeries>>","Cross-correlation diverged,"
    #  ,"may be due to non-converging input (i.e. convex or constant series)."
    #  ,"Series left unshifted."))
    #   delta = 0
    # }
    # else{
    #   delta <- as.vector(cc$lag)[ccmax]
    #   if (cc$acf[ccmax] < threshold) delta = 0
    #   if (shift) return(shiftSeries(y,-delta,padding=padding,using=using,from=from))
    # }
    # return(-delta)
  # }))
}
