#' Optimizes the fit of the given signals to the given series
#' 
#' Calls JS ml-nmr-processing.optimizeSignals from a JS environment built using the V8 package
#' Warning: may crash if given negative ppm. Fix coming.
#' @param x numeric vector, the indepedent variable (ppm)
#' @param y numeric vector, the dependent variable (intensity)
#' @param signals a list or any other object that parses to a valid array of signals for nmr-processing.optimizeSignals on the JS side. 
#' Run about.Signal() for details.
#' @param v8Context the V8 context to be used to run the JS side.
#' You should not need to change this.
#' @param options a list that parses to a valid options argument for 
#' nmr-processing.optimizeSignals on the JS side
#' Run about.signalOptimizationOptions() for details (TBD)
#' @path character indicating where in the context we find ml-nmrProcessing
#' You should not need to change this.
#' @frequency The frequency of the NMR experiment, in Hz
#' @integration character, the method to be used to integrate the signals. See
#' integral
#' @returns a list with "signals" equal the fitted signals, "x" the argument x
#' passed to the function, "fittedY" the total predicted intensity obtained
#' by interpolating the fitted signals on "x" and adding them up, and 
#' "integrals" the integrals of the fitted signals.
#' @import V8
#' @export
optimizeSignals <- function(x,y,signals,v8Context=ct,path="nmrProcessing"
                            ,options=list(optimization=list(kind="lm")
                                          #,linewidth=.1
                                          ,shape=list(kind="pseudoVoigt"
                                                      #,fwhm=.5
                                                   #    ,mu=.7
                                                 ))
                            ,frequency=400,integration="fwhm"){
  if (is.null(options$simulation)){
    options$simulation <- list()
    options$simulation$frequency <- frequency
  } 
  else{
    if (is.null(options$simulation$frequency))
      options$simulation$frequency <- frequency
  }
  
  v8Context$assign("options",options)
  v8Context$assign("signals",signals)
  fittedSignals <- xyApplyJS(x,y,script=paste(path
                                              ,"optimizeSignals(xy,signals,options)"
                                              ,sep="."),v8Context=v8Context)
  fittedXY <- signalsToY(x,fittedSignals,frequency)
  integrals <- unlist(lapply(fittedSignals$peaks, function(signal)
    integral(x,y,signal,method=integration)))
  list(signals=fittedSignals, x=x, fittedY=fittedXY, expY=y, integrals=integrals)
}