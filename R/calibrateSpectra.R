#' Aligns series to a reference
#' Computes the shift that maximizes each series' cross-correlation to the
#' @param Y matrix, series in rows
#' @param ref character, numeric or function. Specifies the reference for
#' alignment. May be a series of the same length as the input series, the index
#' of the row of Y to be used as reference, a function to compute the
#' reference from the input, or a character representing an implemented
#' reference. Valid character references are 'median' (default: reference is
#' the median of the input) and 'mean'. If ref is a function, it run on Y
#' and the result is used as reference.
#' @param shifts logical. If TRUE, returns the shifts; if FALSE, returns the
#' shifted series. Default: FALSE
#' @param ... additional parameters
#' @returns either a numeric vector of shifts or a matrix with the shifted
#' series in its rows
#' @importFrom stats ccf median
alignSeries <- function(Y, ref="median", shifts=FALSE, ...){
  #QC input Y
  if(!is.matrix(Y))
  cat(crayon::red("alignSeries >>",
                  "Y must be a matrix\n"))
  stop()
  #Parse ref and build reference spectrum
  if (is.function(ref)){
    ref <- ref(Y)
  }
  else{
    if (is.character(ref)){
      if (ref=="median") ref <- apply(Y,2,median)
      else if (ref=="mean") ref <- apply(Y,2,mean)
      #TBD?: "centroid" of the set of spectra
      else {
        cat(crayon::red("alignSeries >>",
                        "Invalid character ref: must be 'median; (default) or 'mean'\n"))
        stop()
      }
    }
    else{
      #TBD: numeric reference may be interpreted as index or as a reference series
      if (is.numeric(ref)){
        if (length(ref) == 1) if (ref <= dim(Y)[1]) ref <- Y[ref,]
      }
      else{
        if (length(ref) != dim(Y)[2]) {
          cat(crayon::red("alignSeries >>",
                          "Invalid ref: reference spectrum does not match the length of input spectra\n"))
          stop()
        }
      }
    }
  }
  if (length(ref) != dim(Y)[2]) {
    cat(crayon::red("alignSeries >>",
                    "Invalid ref: filter length does not match the length of input spectra\n"))
    stop()
  }

  t(apply(Y,1,function(y){
    cc <- ccf(y, ref, type="correlation", plot = FALSE, ...)
    ccmax <- which.max(cc$acf)
    shift <- as.vector(cc$lag)[ccmax]
    if (shifts) return(shift)
    trail_l <- if (shift <= 0) 0 else shift
    trail_r <- if (shift <= 0) -shift else 0
    c(rep(0, trail_r), y[(trail_l+1):(length(y)-trail_r)], rep(0,trail_l))
  }))
}

#' Spectra calibration
#' Calibrates spectra by aligning a reference signal
#' @param x numeric, spectra scale e.g. ppm
#' @param Y matrix, intensities, spectra in rows
#' @param matrix character, optional. The sample matrix; either 'urine', 'plasma' or
#' 'other'. Either `matrix` or `RoRef` must be specified.
#' @param RoRef numeric or logical, optional. Specification of the Region of
#' Reference corresponding to the reference signal for calibration. It may be
#' specified as a vector of x coordinates, as a length 2 vector declaring the
#' limits of the region, or as logical filter for x. Either `RoRef` or `matrix`
#' must be provided. If `RoRef` is not provided, the standard reference peak
#' for the given `matrix` will be used as Region of Reference.
#' @param ... additional arguments for `alignSeries` (see details)
#' @details crops the Region of Reference of the spectra matrix and passes it
#' to `alignSeries`, which calculates the shifts that align the spectra on the
#' Region of Reference. Then, it shifts the full spectra by the
#' corresponding amounts.
#' @returns calibrated spectra matrix
#' @importFrom stats ccf
calibrateSpectra <- function(x, Y, RoRef, matrix, ...){
  RoRefs <- list("urine"=c(0.938,0.945),
                 "plasma"=c(),
                 "other"=c(-0.05,0.05))
  #Parse RoRef and construct RoRef filter
  if (missing(RoRef)){
    if (missing(matrix)) {
      cat(crayon::red("alignSeries >>",
                      "Provide 'RoRef' or 'matrix' argument"))
      stop()
    }

    RoRef <- RoRefs[[matrix]]
    RoRef <- x >= RoRef[1] & x <= RoRef[2]
  }
  else{
    if (is.numeric(RoRef)){
      if (length(RoRef==2)){
        RoRef <- x >= RoRef[1] & x <= RoRef[2]
      }
      else{
        RoRef <- x %in% RoRef
      }
    }
    else{
      if (is.logical(RoRef)){
        if (length(x) != length(RoRef))
          cat(crayon::red("alignSeries >>",
                          "RoRef filter does not match x"))
        stop()
      }
      else{
        cat(crayon::red("alignSeries >>",
                        "Invalid roRef: must be logical, numeric or 'matrix'"))
        stop()
      }
    }
  }
  #Align on RoRef and get shifts
  shifts <- alignSeries(Y[,RoRef], shifts=TRUE, ... )
  #Shift whole spectra by the corresponding shifts
  t(sapply(1:dim(Y)[1],function(i){
    shift <- shifts[i]
    y <- Y[i,]
    trail_l <- if (shift <= 0) 0 else shift
    trail_r <- if (shift <= 0) -shift else 0
    c(rep(0, trail_r), y[(trail_l+1):(length(y)-trail_r)], rep(0,trail_l))
  }))
}
