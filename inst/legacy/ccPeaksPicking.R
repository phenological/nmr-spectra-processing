#peaks is a list or data.frame (for compatibility with autoPeaksPicking)
#list is coerced to data.frame using peaks %>% toJSON %>% fromJSON
#each peak is specified by x,y,shape = {kind,fwhm,mu}
#x is the mean of the peak, y is its maximum; labels x y are used for compatibility
#with autoPeaksPicking
#only x is mandatory; defaults: y=1 shape={kind="pseudoVoigt",fwhm=0.5/frequency,mu=0.7}
ccPeakPicking <- function(x,y,peaks,frequency){
  #parse and qc input
  if (length(x)!=length(y)){
    cat(crayon::red("nmrSpectraProcessing::ccPeakPicking >>"
    ,"x and y must be of the same length"
    ))
    stop()
  }
  if (is.list(peaks))
    peaks <- jsonlite::fromJSON(jsonlite::toJSON(peaks,auto_unbox=TRUE))
  for (i in 1:dim(peaks)[1]){
    
  }
  
}