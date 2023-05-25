#' Interpolates fitted signals on the given set of points and combines them to 
#' generate the overall intensity
#' 
#' @param x numeric, points to interpolate
#' @param signals data.frame with columns x, y, width, shape[kind,fwhm], ... 
#' A ml-nmrProcessing array of signals parsed into R through a V8 context. 
#' Run about.signal() for details.
#' @param frequency numeric, the frequency in Hz of the NMR experiment
#' @returns a numeric vector of intensities computed on x
#' @export
signalsToY <- function(x,signals,frequency=400){
  res <- 0
  #trick to make sure signals is a data.frame even if signalsToY received a list
  if(!is.data.frame(signals) & is.list(signals)){
    signals <- jsonlite::fromJSON(jsonlite::toJSON(signals))
  }
  for (signal in signals$peaks){
    for (i in 1:dim(signal)[1]){
      peaks <- signal[i,]
      if (peaks$shape$kind[1] == "gaussian")
        res <- res + rowSums(gaussian(x,max=peaks$y,mean=peaks$x
                                      ,fwhm=peaks$shape$fwhm/frequency))
      else{
        if(peaks$shape$kind[1] == "lorentzian"){
          res <- res + rowSums(lorentzian(x,max=peaks$y,mean=peaks$x
                                          ,fwhm=peaks$shape$fwhm/frequency))
        }
        else{
          if(peaks$shape$kind[1] == "pseudoVoigt")
            res <- res + rowSums(pseudoVoigt(x,max=peaks$y,mean=peaks$x
                                             ,fwhm=peaks$shape$fwhm/frequency
                                             ,mu=peaks$shape$mu))
          else{
            cat(crayon::red("signalsToY>>","unknown kind in signal",i,"\n"))
          }
        }
      }
    }
  }
  return(res)
}