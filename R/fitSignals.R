#' Picks peaks ands fits signals
#' 
#' @export
#Built flexible so that we can try different peak picking and
#optimization algorithms, but at the moment only uses autoPeaksPicking and
#optimizeSignal.
fitSignals <- function(x,y,expectedPeaks,roi,signals,grouping
                       ,v8Context=ct,path="nmrProcessing",frequency=400
                       ,integration="fwhm"
                       ,peakPicking=autoPeaksPicking
                       ,peakPickingOptions=list(shape=list(kind="pseudoVoigt"
                                                           ,fwhm=0.5
                                                           ,mu=0))
                       ,optimize=TRUE,optimization=optimizeSignals
                       ,optimizationOptions=list(optimization=list(kind="lm")
                                                 ,shape=list(kind="pseudoVoigt"
                                                             ,fwhm=0.5
                                                             ,mu=0))){
  #Quality check on x y length
  if (length(x) != length(y)){
    cat(crayon::red("fitSignals>>","x and y length do not match\n"))
    stop()
  }
  #Parse and apply the roi filter
  if (!missing(roi)){
    if (!is.logical(roi)){
      if (is.numeric(roi)) roi <- x>=roi[1] & x<=roi[2]
      else
        stop("roi must be either logical or numeric (extremes of the region of interest)")
    }
    y <- y[roi]
    x <- x[roi]
  }
  if (is.null(optimizationOptions$simulation))
    optimizationOptions$simulation = list(frequency=frequency)
  else{
    if (is.null(optimizationOptions$simulation$frequency))
      optimizationOptions$simulation$frequency = frequency
  }
  if (is.null(peakPickingOptions$frequency))
    peakPickingOptions$frequency=frequency
  
  #Pick peaks
  peaksPicked <- autoPeaksPicking(x,y,v8Context=v8Context,frequency=frequency
                                  ,path=path,options=peakPickingOptions)
  #Check and adjust results of peak picking
  if (is.null(peaksPicked)){
    cat(crayon:yellow("fitSignals>>","No peaks found\n"))
    return(NULL)
  }
  if (countPeaks(peaksPicked) < expectedPeaks){
    cat(crayon::yellow("fitSignals>>","Not enough peaks found\n"
                       , "Returning peaks picked\n"))
    return(peaksPicked)
  }
  if(countPeaks(peaksPicked) > expectedPeaks){
    cat(crayon::yellow("Too many peaks found. Dropping lowest intensity peaks\n"
                       ,"This is unlikely to solve your issues\n"
                       ,"You may need to modify your model of the region\n"))
    while(countPeaks(peaksPicked) > expectedPeaks){
      bye <- which.min(peaksPicked$y)
      peaksPicked <- peaksPicked[-bye,]
    }
  }
  #Check and apply peaks-to-signals grouping
  if (missing(grouping)){
    #If neither grouping nor signals are provided, assigns all picked peaks
    #to a single signal
    if(missing(signals)){
      signals <- peaksToSignals(peaksPicked,list(1:(dim(peaksPicked)[1])))
    }
    else{
      #If grouping is not provided and a single signal is provided,
      #all picked peaks are assigned to that signal. Signal's $peaks
      #are overwritten if previously existing. Delta assigned to mean of
      #picked peaks if not given.
      if(length(signals==1)){
        if (is.null(signals[[1]]$delta)) signals[[1]]$delta <- mean(peaksPicked$x)
        signals[[1]]$peaks <- peaksPicked
      }
      else
        stop("peaks grouping into signals must be specified")
    }
  }
  else{
    if (length(grouping)!=length(signals))
      stop("Number of signals does not match number of peak groups")
    else{
      #assigned peaks to signals according to the given grouping. Overwrite
      #peaks if already existing in the signal. Set delta to mean of grouped
      #peaks if not given in the signal.
      signalsFromPeaks <- peaksToSignals(peaksPicked,grouping)
      for (i in 1:length(signals)){
        if (is.null(signals[[i]]$delta)) signals[[i]]$delta <- signalsFromPeaks[i,]$delta
        signals[[i]]$peaks <- signalsFromPeaks[i,]$peaks
      }
    }
  }
  #Optimize signals and format output as required
  if (optimize)
    optimization(x, y, signals, v8Context=ct, path="nmrProcessing"
                 , options=optimizationOptions, frequency=frequency
                 , integration=integration)
  else{
    fittedY <- signalsToY(x,signals,frequency)
    integrals <- unlist(lapply(signals$peaks, function(signal)
      integral(x,y,signal,method=integration)))
    list(signals=signals, x=x, fittedY=fittedY, expY=y, integrals=integrals)
  }
}