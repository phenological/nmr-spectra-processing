#' Automatic peak picking on a region of interest using Covas' method
#' 
#' @export
#' Peak picking runs on the JS side, using nmrProcessing.autoPeaksPicking.
#' Communication between R and JS is mediated by jsonlite. 
#' See xyApplyJS for details.
#' @param x numeric, x coordinate values (ppm)
#' @param y numeric, y coordinate values (intensity)
#' @param v8Context the V8 context used to communicate with JS. You should not
#' need to change this. See xyApplyJS for details
#' @param path character, module of the JS bundle where nmrProcessing is located.
#' You should not need to change this.
#' @options Options for nmrProcessing.autoPeaksPicking. 
#' Run about.autoPeaksPicking() for details (TBD)
#' @returns a data.frame of peaks. Run about.NMRPeak1D() for details.
#' @import V8
autoPeaksPicking <- function(x,y,v8Context=ct,path="nmrProcessing",frequency=400
                             ,options=list(shape=list(kind="pseudoVoigt"
                                                      ,fwhm=0.5
                                                      ,mu=0))){
  if (is.null(options$frequency)) options$frequency <- frequency
  v8Context$assign("options",options)
  xyApplyJS(x,y,script=paste(path,"xyAutoPeaksPicking(xy,options)",sep="."))
}