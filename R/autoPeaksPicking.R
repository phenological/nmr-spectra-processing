#' Automatic peak picking on a region of interest using Covas' method
#' 
#' Peak picking runs on the JS side, using nmrProcessing.autoPeaksPicking.
#' Communication between R and JS is mediated by jsonlite. 
#' See xyApplyJS for details.
#' @param x numeric, x coordinate values (ppm)
#' @param y numeric, y coordinate values (intensity)
#' @param v8Context the V8 context used to communicate with JS. You should not need to change this. See xyApplyJS for details
#' @param path character, module of the JS bundle where nmrProcessing is located. You should not need to change this.
#' @options Options for nmrProcessing.autoPeaksPicking. Run about.autoPeaksPicking() for details (TBD)
#' @returns a data.frame of peaks. Run about.NMRPeak1D() for details.
#' @import V8
#' @export
autoPeaksPicking <- function(x,y,v8Context=ct,path="nmrProcessing"
                             ,minimum = 1
                             ,recursive=FALSE
                             ,threshold=3,rEnd=1,rStep=0.5
                             ,frequency=400
                             ,options=list(shape=list(kind="pseudoVoigt"
                                               #       ,fwhm=0.5
                                               #       ,mu=.7
                                               ))
                             ,incomplete=FALSE){
  if (is.null(options$frequency)) options$frequency <- frequency
  if (is.null(options$thresholdFactor)) options$thresholdFactor <- threshold
  v8Context$assign("options",options)
  pickedPeaks <- xyApplyJS(x,y,script=paste(path,"xyAutoPeaksPicking(xy,options)"
                                            ,sep="."))
  cp <- countPeaks(pickedPeaks)
  if (cp < minimum){
    if (recursive){
      if (threshold < rEnd){
        crayon::yellow("autoPeaksPicking>>"
                       ,"Found",cp
                       ,"peaks, less than the specified minimum of ", minimum)
        if (incomplete & cp > 0) return(pickedPeaks) 
        return(NULL)
      }
      #cat(threshold,options$thresholdFactor)
      tf <- options$thresholdFactor - rStep
      options$thresholdFactor <- tf
      #cat(tf,options$thresholdFactor,"\n")
      # print("recurring")
      return(autoPeaksPicking(x,y,v8Context,path=path
                              ,minimum=minimum,recursive=TRUE
                              ,threshold=tf,rEnd=rEnd,rStep=rStep
                              ,frequency=frequency,options=options
                              )
             )
    }
    else{
      if (incomplete & (cp > 0)) return(pickedPeaks)
      return(NULL)
    } 
  }
  else return(pickedPeaks)
}