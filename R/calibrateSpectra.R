#Some utils for working with fusion::NMRSignal1D
normalizeSignal <- function(signal){
  summit <- max(sapply(signal@peaks,function(aPeak) aPeak@y))
  signal@peaks <- lapply(signal@peaks, function(aPeak){
    aPeak@y <- aPeak@y / summit
    return(aPeak)
  })
  return(signal)
}

#Signal domain is defined as widthFactor around the external peaks
#widthFactor=3 makes sense in concept, but very larger domains are necessary for
#applications such as calibration
signalDomain <- function(aSignal, widthFactor = 3){
  offset <- aSignal@shape$params$fwhm * widthFactor * c(-1,1)
  peaksRange <- range(sapply(aSignal@peaks, function(aPeak) aPeak@x))
  peaksRange + offset
}

#Auxiliary function, returns the shift values but does not shift so that it can
#be used by the actual exports to either calibrate the spectra to a signal or
#calibrate a signal to a spectrum
calibrateToSignal <- function(x, Y, signal, rOref=signalDomain(signal,30)
                              , maxShift=1/3,threshold=0.2, ...){
  rOref <- crop(x,roi=rOref)
  x <- x[rOref]
  #Align normalized spectra on rOref to reference and get shifts
  if (is.matrix(Y)){
    normS <- t(apply(Y[,rOref],1,function(y) y / max(y)))
  } else{
    normS <- Y[rOref]
    normS <- normS / max(normS)
  }
  alignSeries(normS, signalToY(normalizeSignal(signal),x),shift=FALSE
              ,lag.max=length(x) * maxShift, threshold=threshold, ...)
}

#' Signal calibration
#' Calibrates a (model) signal to the experimental spectrum
#' @param x numeric, spectra scale e.g. ppm
#' @param y numeric, spectrum intensities
#' @param signal NMRSignal1D
#' @param rOref numeric, optional, limits of the Region of Reference within which the reference signal will be calibrated. The default encompasses a region 30 linewidths around the external peaks.
#' @param maxShift numeric, the maximum allowed shift as a fraction of the rOref, see ccf, lag.max
#' @param threshold numeric, minimum cross-correlation requiered for alignment. See alignSeries
#' @param ... additional arguments for alignSeries, see Details
#' @returns a NMRSignal1D with peak maxima shifted to maximize cross-correlation between its trace and the reference spectrum
calibrateSignal <- function(x,y,signal, maxShift=1/3,threshold=0.2
                            ,rOref=signalDomain(signal,30), ...){
  aShift <- calibrateToSignal(x,y,signal,rOref, maxShift=1/3
                              ,threshold=0.2,...)
  aShift <- aShift * (x[2] - x[1])
  for (i in 1:length(signal@peaks)){
    signal@peaks[[i]]@x <- signal@peaks[[i]]@x - aShift
  }
  signal@chemicalShift <- signal@chemicalShift - aShift
  return(signal)
}

#' Spectra calibration
#' Calibrates spectra
#' @param x numeric, spectra scale e.g. ppm
#' @param Y matrix, intensities, spectra in rows
#' @param ref numeric or NMRSignal1D, the reference signal for calibration. Either 'glucose', 'alanine', 'tsp', or a fusion::NMRSignal1D
#' @param frequency numeric, the spectrometer frequency in MHz, only important for character ref. Default: 600
#' @param rOref numeric, optional. Limits of the Region of Reference within which the reference signal for calibration will be aligned. Use if the calibrated chemical shift is too far from the experimental signal, or if you need to avoid interfering signals.
#' @param j numeric, optional. Experimental coupling constant. Only relevant for ref="alanine" or "glucose". If provided, it replaces the internal model's coupling constant. Can be used to improve calibration to alanine and glucose, which may not be accurate due to fluctuations in the observed j.
#' @param maxShift numeric, the maximum allowed shift as a fraction of the rOref, see ccf, lag.max
#' @param threshold numeric, minimum cross-correlation requiered for alignment. See alignSeries
#' @param padding character, the method to be used to fill the empty extremes of the shifted spectra. See 'pad' for details. Default: "sampling".
#' @param from numeric, optional. Filter for the spectral region to be used for padding after the spectrum is shifted in the "sampling" method. See 'pad' for details. Default: last 1/15th points
#' @param plot logical. Argument passed to ccf through alignSeries; if TRUE, each spectrum's correlation to the reference is plotted. Default: FALSE. You should not need to change the default unless you are getting misalignments and you want to check the cross-correlation to fine-tune the alignment parameters.
#' @param ... additional arguments for alignSeries, see Details
#' @details crops the Region of Reference of the spectra matrix and passes it to alignSeries, along with a model of the reference signal, to calculate the shifts that align the spectra on the rOref to the ref signal. Then, it shifts the full spectra by the corresponding amounts.
#' @returns calibrated spectra matrix
#' @importFrom stats ccf
calibrateSpectra <- function(x, Y,ref=c("tsp","glucose","alanine"
                                        ,"NMRSignal1D see documentation")[1]
                             ,frequency=600, maxShift=1/3
                             ,threshold=0.2, padding="zeroes"
                             ,from=as.integer(length(x)*14/15):length(x)
                             ,rOref, j, ...){
  
  if (missing(j)) j <- NULL else j <- j/2
  
  make.ref <- function(ref, frequency, j){
    if (ref == "tsp"){
      return(new("NMRSignal1D"
                 ,id = "TSPbration"
                 ,chemicalShift = 0
                 ,shape=list(name = "pseudoVoigt"
                             ,params=list(mu=.85,fwhm=0.0001,base=0))
                 ,peaks = list(
                   new("NMRPeak1D", x=0, y=1)
                 )
      ))
    }
    if (ref == "alanine"){
      if (is.null(j)) j <- 3.63/frequency
      return(new("NMRSignal1D"
                 ,id = "Alabration"
                 ,chemicalShift = 1.48
                 ,shape=list(name = "pseudoVoigt"
                             , params=list(mu=.85,fwhm=.6/frequency,base=0))
                 ,peaks = list(
                   new("NMRPeak1D", x=1.48 - j, y=1)
                   ,new("NMRPeak1D", x=1.48 + j, y=1)
                 )
      ))
    }
    if (ref == "glucose"){
      if (is.null(j)) j <- 1.866 / frequency
      return(new("NMRSignal1D"
                 ,id = "Glucobration"
                 ,chemicalShift = 5.223
                 ,shape=list(name = "pseudoVoigt"
                             , params=list(mu=.85,fwhm=.6/frequency,base=0))
                 ,peaks = list(
                   new("NMRPeak1D", x=5.223 - j, y=1)
                   ,new("NMRPeak1D", x=5.223 + j, y=1)
                 )
      ))
    }
    cat(crayon::red("nmr.spectra.processing::calibrateSpectra >>"
                    ,"invalid reference\n"))
    stop()
  }
  
  #create ref if needed
  if (is.character(ref)) ref <- make.ref(ref, frequency, j)
  if (!("NMRSignal1D" %in% is(ref))){
    cat(crayon::red("nmr.spectra.processing::calibrateSpectra >>"
                    ,"invalid reference\n"))
    stop()
  }
  
  #compute rOref if needed (tsp is handled differently because it's modeled as a "line")
  if (missing(rOref)){
    if (ref@id == "TSPbration"){
      rOref <- c(-0.02,0.02)
    } else rOref <- signalDomain(ref)
  }
  
  #Align normalized spectra on rOref to reference and get shifts
  shifts <- calibrateToSignal(x,Y,ref,rOref, maxShift=1/3
                              ,threshold=0.2,...)
  
  #Shift whole spectra by the corresponding shifts
  if (is.matrix(Y)){
    t(sapply(1:dim(Y)[1],function(i){
      shift <- shifts[i]
      y <- Y[i,]
      shiftSeries(y,shift, padding=padding, from=from)
    }))
  } else{
    shiftSeries(y,shift,padding=padding,from=from)
  }
}
