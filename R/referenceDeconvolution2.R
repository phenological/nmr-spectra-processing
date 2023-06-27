#Candidate to utilities.R
crop <- function(x, roi){
  x >= roi[1] & x <= roi[2]
}

targetFunction <- function(parS, xx) {
  mu = parS[[1]]
  fwhm = parS[[2]]
  base = parS[[3]]
  nSignals = (length(parS) - 3)/ 2
  signal = base
  for (i in 1:nSignals) {
    signal = signal + pseudoVoigt(xx, mu=mu, mean=parS[[i* 2 + 2]], max=parS[[i * 2 + 3]], fwhm=fwhm)
  }
  return(signal)
}

#' Reference deconvolution using TSP peak as reference
#' @importFrom gsignal hilbert
#' @export
#' 
referenceDeconvolution2 <- function(ppm,
                                    y,
                                    rOref=c(-0.02,0.02),
                                    nmrSignal1D = NULL,
                                    padding="zeroes",
                                    using=c(9.5,10),
                                    from=y,
                                    noise=c(0, 0.5)) {
  #inverse fourier transform with normalization
  ifft <- function(x){
    fft(x, inverse = TRUE) / length(x)
  }
  
  from <- rOref[[1]]
  to <- rOref[[2]]
  
  if(is.null(nmrSignal1D)) {
    shape = list("name" = "pseudoVoigt", "params" = list("mu"=0.85, 
                                                         "fwhm"=0.00105, 
                                                         "minmu"=0, 
                                                         "maxmu"=1,
                                                         "minfwhm"=0,
                                                         "maxfwhm"=0.00105 + 0.003,
                                                         "base" = 0))
    
    nmrSignal1D <- new("NMRSignal1D", peaks=list(new("NMRPeak1D", x = -0.0054, y = 1, fwhm = 0.002), 
                                                 new("NMRPeak1D", x = 0, y = 100, fwhm = 0.002), 
                                                 new("NMRPeak1D", x = 0.0054, y = 1, fwhm = 0.002)), 
                       shiftRange = 0.0005, 
                       shape = shape)
    
  }
  
  fi <- ppm < to & ppm > from
  ppm_ = ppm[fi]
  model = new("NMRSignalModel", signalsInput = list(nmrSignal1D),
              from = from,
              to =  to,
              ppm = ppm_,
              experimental = y[fi])
  
  model = fitModel(model)
  fitted = model@fitted
  signalOutput = model@signalsOutput[[1]]
  
  maxSatelite = (signalOutput@peaks[[1]]@y +
                   signalOutput@peaks[[3]]@y) / 2
  
  signalOutput@peaks[[1]]@y <- maxSatelite
  signalOutput@peaks[[3]]@y <- maxSatelite
  mu = signalOutput@shape$params$mu
  fwhm = signalOutput@shape$params$fwhm
  
  fitted = 0
  for (peak in signalOutput@peaks)  {
    fitted = fitted + pseudoVoigt(ppm_, mu=mu, mean=peak@x, max=peak@y, fwhm=fwhm)
  }
  
  
  
  #noise = runif(length(ppm), noise[[1]], noise[[2]]) * 0
  refModel = rep(0, length(ppm))
  refModel[fi] = fitted #+ noise[fi]
  ref = rep(0, length(ppm))
  ref[fi] = y[fi] #+ noise[fi]
  for( i in which(fi)) {
    ref[[i]] =  (y[i - 1] +  y[i] +  y[i + 1])/3
  }
  #move to time domain
  ref <- ifft(hilbert(ref))
  refModel <- ifft(hilbert(refModel))
  y <- ifft(hilbert(y))
  #deconv., correct and go back to freq. domain
  return (Re(fft(y * refModel / ref)))
}

# for (index in 1:20 ) {
#   rOref <- c(-0.01,0.01)
#   
#   fi <- ppm > rOref[1] & ppm < rOref[2]
#   y2 <- referenceDeconvolution2(ppm, Y[index,], rOref=rOref, noise=c(0, 0.1))
#   plot(ppm[fi], y2[fi], type = 'l')
#   lines(ppm[fi], Y[index, fi], type = 'l', col='red')
#   
#   fi <- ppm > 11.9 & ppm < 12.1
#   plot(ppm[fi], y2[fi], type = 'l')
#   lines(ppm[fi], Y[index, fi], type = 'l', col='red')
#   
#   fi <- ppm > 3.45 & ppm < 3.55
#   plot(ppm[fi], y2[fi], type = 'l')
#   lines(ppm[fi], Y[index, fi], type = 'l', col='red')
#   
#   fi <- ppm > 8.15 & ppm < 8.25
#   plot(ppm[fi], y2[fi], type = 'l')
#   lines(ppm[fi], Y[index, fi], type = 'l', col='red')
# }

