#' Constructs a list of signals from a data.frame of peaks
#' 
#' Useful for optimizing the results of autoPeaksPicking
#' @param peaks data.frame with peaks to be assigned to signals. Run about.NMRPeak1D() for details.
#' @param grouping list. Specification of how the peaks are to be partitioned in signals. Each element of the list is a numeric vector giving the row indices of the peaks that will conform each signal. Run about.Signal() for details.
#' @returns a list of signals with "delta" equal the mean x of assigned peaks, and "peaks" equal the data.frame of assigned peaks.
#' @import V8
#' @export
#TBD: implement optional arguments to pass more attributes to the signals
#(not prioritary, maybe even a bad idea). Allow grouping to be specified through
#peak ids (have not met a scenario where it would be needed but sounds plausible)
peaksToSignals <- function(peaks,grouping){
  #trick to make sure that peaks is a data.frame even if signalsToY received a list
  #Obsolete and dangerous?
  if(!is.data.frame(peaks) & is.list(peaks)){
    peaks <- jsonlite::fromJSON(jsonlite::toJSON(peaks))
  }
  
  cp <- countPeaks(peaks)
  if(cp == 0){
    cat(crayon::yellow("peaksToSignals>>", "no peaks to process"))
    return(NULL)
  }
  #If no grouping, all peaks to one signal (singulet)
  if(missing(grouping))
    grouping <- list(1:cp)
  qc <- sort(unlist(grouping))
  if (length(qc)!=cp | sum(qc==(1:length(qc))) != length(qc)){
    cat(crayon::red("peaksToSignals>>"
                    ,"provided",countPeaks(peaks),"peaks, but groups contain"
                    ,length(qc), "peaks\n"
                    ,"All expectedPeaks must be grouped into the specified signals\n"))
    stop()
  }
  #lapply(grouping, function(grp){
  #  delta <- mean(peaks[grp,]$x)
  #  pks <- peaks[grp,]
  # list(delta=delta,peaks=pks) 
  #})
  #jsonlite::fromJSON(jsonlite::toJSON(signals,auto_unbox=TRUE))
  delta <- sapply(grouping, function(group){
    mean(peaks[group,]$x)
  })
  picos <- lapply(grouping, function(group){
    peaks[group,]
  })
  df <- data.frame(delta=delta)
  df$peaks <- picos
  df
}