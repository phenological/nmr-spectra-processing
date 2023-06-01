#' Picks peaks ands fits signals
#' 
#' @import V8
#' @export
#Built flexible so that we can try different peak picking and
#optimization algorithms, but at the moment only uses autoPeaksPicking and
#optimizeSignal.
fitSignals <- function(x,y,roi,signals,grouping,dropExcess=.75
                       ,v8Context=ct,path="nmrProcessing",frequency=400
                       ,integration="fwhm"
                       ,peakPicking=autoPeaksPicking
                       ,recursivePicking=TRUE
                       ,threshold=3
                       ,rEnd=1
                       ,rStep=0.5
                       ,peakPickingOptions=list(shape=list(kind="pseudoVoigt"
                                                       #    ,fwhm=0.5
                                                           ,mu=.7))
                       ,optimize=TRUE,optimization=optimizeSignals
                       ,optimizationOptions=list(optimization=list(kind="lm")
                                                 ,linewidth=0.1
                                                 ,shape=list(kind="pseudoVoigt"
                                                          #   ,fwhm=0.5
                                                             ,mu=.7))){
  #Quality check on x y length
  if (length(x) != length(y)){
    cat(crayon::red("fitSignals>>","x and y length do not match\n"))
    stop()
  }
  #Parse and apply the roi filter
  if (!missing(roi)){
    if (!is.logical(roi)){
      if (is.numeric(roi)){
        #Test: if i let xyAutoPeaskPicking.ts look at the whole spectrum
        #before zooming into the target zone, maybe it will do a better
        #work at estimating noise
        #peakPickingOptions$from = roi[1]
        #peakPickingOptions$to = roi[2]
        roi <- x>=roi[1] & x<=roi[2]
      } 
      else
        stop("roi must be either logical or numeric (extremes of the region of interest)")
    }
    y <- y[roi]
    x <- x[roi]
  }
  #Add frequency parameter to options as needed
  if (is.null(optimizationOptions$simulation))
    optimizationOptions$simulation = list(frequency=frequency)
  else{
    if (is.null(optimizationOptions$simulation$frequency))
      optimizationOptions$simulation$frequency = frequency
  }
  if (is.null(peakPickingOptions$frequency))
    peakPickingOptions$frequency=frequency
  
  #compute minimum expected number of peaks, if possible
  if (missing(grouping)){
    #no grouping, no signals or only 1 signal -> anything is ok
    if (missing(signals)) expectedPeaks <- 1
    else{
      if (length(signals) == 1) expectedPeaks <- 1
      else{
        #no grouping, multiple signals -> error
        cat(crayon::red("fitSignals>>","cannot group peaks into multiple signals",
                        "without grouping"))
        stop()
      }
    }
  }
  else{
    #grouping <- expect peaks specified in grouping
    expectedPeaks <- length(unlist(grouping))
  }
  
  #Pick peaks
  peaksPicked <- autoPeaksPicking(x,y,minimum=expectedPeaks
                                  ,v8Context=v8Context,frequency=frequency
                                  ,path=path
                                  ,recursive=recursivePicking
                                  ,threshold=3 ,rEnd=1 ,rStep=0.5
                                  ,options=peakPickingOptions)
  cp <- countPeaks(peaksPicked)
  
  ##No peaks picked --> return NULL
  if (cp==0){
    cat(crayon::yellow("fitSignals>>","Not enough peaks picked\n"))
    return(NULL)
  }
  
  ##Invalid attribute names in peaks picked --> return NULL
  if (sum(is.na(peaksPicked$x))>0 | sum(is.na(peaksPicked$y))>0) {
    cat(crayon::yellow("fitSignals>>","Got NA in peak attributes!!!!\n"))
    return(NULL)
  }
  
  #no grouping -> update expectedPeaks to the number picked
  if (missing(grouping)){
    expectedPeaks <- cp
    }
          
  ##More peaks picked than expected -> drop until OK, if requested on call
  ##Otherwise return NULL
  dropped <- data.frame()
  if (cp > expectedPeaks){
    cat(crayon::yellow("fitSignals>>", cp
                       , "peaks found; that is too many\n"))
    if (((cp - expectedPeaks) / expectedPeaks) <= dropExcess){
      cat(crayon::yellow("fitSignals>>"
                         ,"Dropping lowest intensity peaks\n"
                         ,"This is unlikely to solve your issues\n"
                         ,"You may need to modify your model of the region\n"
                         ))
      while(cp > expectedPeaks){
        bye <- which.min(peaksPicked$y)
        dropped <- rbind(dropped,peaksPicked[bye,])
        peaksPicked <- peaksPicked[-bye,]
        cp <- cp - 1 #countPeaks(peaksPicked)
      }
    }
    else return(NULL)
  }
  
  #Fill in signals and grouping as required and update signals
  #with grouped peaks' data
  #Defaults were not defined on arguments as they depend on the number of
  #picked Peaks
  if (missing(signals)){
    if (missing(grouping))
      #Neither signals nor grouping -> assume one signal (singlet) per pickedPeak
      signals <- peaksToSignals(peaksPicked,1:cp)
    else
      #grouping without signals -> use grouping to create fresh signals
      #(all singles)
      signals <- peaksToSignals(peaksPicked,grouping)
  }
  else{
    if (missing(grouping)){
      #signals without grouping -> it is already constrained to be
      #a single signal, so we can update it with picked peaks and continue
      signalsFromPeaks <- peaksToSignals(peaksPicked)
      if (is.null(signals[[1]]$delta))
        signals[[1]]$delta <- signalsFromPeaks[1,]$delta
      signals[[1]]$peaks <- signalsFromPeaks[1,]$peaks
    }
    else{
      #both signals an grouping are given -> update signals with peak data
      signalsFromPeaks <- peaksToSignals(peaksPicked,grouping)
      for (i in 1:length(signals)){
        if (is.null(signals[[i]]$delta))
          signals[[i]]$delta <- signalsFromPeaks[i,]$delta
        signals[[i]]$peaks <- signalsFromPeaks[i,]$peaks
      }
    }
  }
  
  #Optimize signals and format output as required
  if (optimize){
    res <- optimization(x, y, signals, v8Context=ct, path="nmrProcessing"
                 , options=optimizationOptions, frequency=frequency
                 , integration=integration)
    res$dropped <- dropped
    res
  }
  else{
    fittedY <- signalsToY(x,signals,frequency)
    integrals <- unlist(lapply(signals$peaks, function(signal)
      integral(x,y,signal,method=integration)))
    list(signals=signals, x=x, fittedY=fittedY
         , expY=y, integrals=integrals, dropped=dropped)
  }
}