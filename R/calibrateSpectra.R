#' Signal calibration
#' 
#' Calibrates a \code{\linkS4class{NMRSignal1D}} to the experimental spectrum
#' @param ppm numeric, chemical shift scale of the spectrum
#' @param y numeric, spectrum intensities.
#' @param signal \code{\linkS4class{NMRSignal1D}} to be calibrated
#' @param rOref numeric, optional, limits of the Region of Reference within
#'  which the reference signal will be calibrated. The default encompasses a
#'  region 30 linewidths around the external peaks.
#' @param maxShift numeric, the maximum allowed shift as a fraction of the rOref,
#'  see \code{\link[stats]{ccf}}.
#' @param threshold numeric, minimum cross-correlation required for alignment.
#'  See \code{\link{alignSeries}}.
#' @param lambda, numeric. Parameter to \code{\link[ptw]{asysm}}. \code{calibrateSignal} uses
#'\code{\link{baselineCorrection}} with an unusually small \code{lambda} to
#' flatten broad peaks that interfere with signal alignment. The default
#'  \code{lambda=1e3} works well in typical cases.
#' @param ... additional arguments for \code{\link{alignSeries}}.
#' @details Interpolates the signal with \code{\link{signalToY}} and aligns the
#' resulting trace to the spectrum in the rOref using \code{\link{alignSeries}}.
#' This is useful e.g. for fine-tuning a quantification model.
#' @returns a \code{\linkS4class{NMRSignal1D}} with peak maxima shifted to
#' align with the reference spectrum
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
#' @importFrom stats ccf
#' @export
calibrateSignal <- function(ppm,y,signal
                            , maxShift=1/3,threshold=0.2
                            ,rOref=signalDomain(signal,30)
                            ,lambda=1e3, ...){
  aShift <- calibrateToSignal(ppm,y,signal,rOref, maxShift=maxShift
                              ,threshold,lambda=lambda,...)
  aShift <- aShift * (ppm[2] - ppm[1])
  for (i in 1:length(signal@peaks)){
    signal@peaks[[i]]@x <- signal@peaks[[i]]@x - aShift
  }
  signal@chemicalShift <- signal@chemicalShift - aShift
  return(signal)
}

#' Spectra calibration
#' 
#' Calibrates spectra
#' @param ppm numeric, spectra chemical shift scale
#' @param Y matrix, intensities, spectra in rows
#' @param ref character or \code{\linkS4class{NMRSignal1D}}, the reference signal 
#' for calibration. Either "glucose", "alanine", "tsp", "serum", or a
#' \code{\linkS4class{NMRSignal1D}}. "tsp" calibrates to a singlet at 0 ppm. "alanine"
#' calibrates to a doublet at 1.48 ppm with j=7.26. "glucose" aligns to a doublet
#' at 5.223 with j=3.63. "serum" aligns to the glucose double with a (typically)
#' extended range of shifting.
#' @param frequency numeric, the spectrometer frequency in MHz. Only used for
#' character \code{ref}. Default: 600
#' @param rOref numeric, optional. Limits of the Region of Reference within
#' which the reference signal for calibration will be aligned. Use if the
#' calibration signal is very shifted, or if you need to avoid interfering signals.
#' Otherwise, the default should be fine.
#' @param cshift, numeric, optional. Chemical shift of the reference. If 
#' provided, it replaces the internal model's chemical shift. Use if you want 
#' to position your signal differently from the standard model (see \code{ref})
#' @param j numeric, optional. Coupling constant of the reference. If provided, 
#' it  replaces the internal model's coupling constant. Use it to fine tune
#' calibration to alanine and glucose, which may not be accurate due to
#'  fluctuations in the observed \emph{j}.
#' @param maxShift numeric, the maximum allowed shift as a fraction of the rOref,
#'  see \code{\link[stats]{ccf}}
#' @param threshold numeric, minimum cross-correlation required for alignment.
#'  See \code{\link{alignSeries}}
#' @param padding character, the method to be used to fill the empty extremes of
#'  the shifted spectra. See \code{\link{pad}} for details. Default: "zeroes".
#' @param from numeric, optional. Filter for the spectral region to be used for
#'  padding after the spectrum is shifted in the "sampling" method. See 
#'  \code{\link{pad}} for details. Default: last 1/15th points
#' @param lambda numeric. You should not need to tamper with this parameter
#' unless you are attempting something extraordinary. See \code{\link{calibrateSignal}}
#' for details
#' @param ... additional arguments for \code{\link{alignSeries}}, see \strong{Details}
#' @details Interpolates the reference signal with \code{\link{signalToY}} to 
#' obtain a trace. Then, it aligns each row of the spectra matrix to this trace,
#'  within the \code{rOref},  using \code{\link{alignSeries}}. Passing a NMRSignal1D
#'  to \code{ref} allows to calibrate to any custom signal; but if you are aligning
#'  to a custom doublet or singlet, it's more practical to set \code{ref} to either
#'  "tsp" (for singlets) or "alanine" (for doublets) and adjust the signal parameters
#'  through the arguments \code{cshift} and \code{j}. This way \code{calibrateSpectra}
#'  handles the annoying task of translating the signal parameters into a NMRSignal1D
#'  object.
#' @returns calibrated spectra matrix
#' @importFrom stats ccf
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
#' @export
calibrateSpectra <- function(ppm, Y,ref=c("tsp","glucose","alanine","serum"
                                        ,"an NMRSignal1D, see documentation")[1]
                             ,frequency=600, maxShift=1/3,threshold=0.2
                             ,rOref, cshift, j
                             ,padding=c("zeroes","circular","sampling")[1]
                             ,from=as.integer(length(ppm)*14/15):length(ppm)
                             ,lambda=1e3 , ...){
  rowLabels <- rownames(Y)
  #Type check and casting
  if (!is.numeric(ppm)){
    cat(crayon::yellow("nmr.spectra.processing::calibrateSpectra >>"
                       ,"Non-numeric argument ppm being cast as.numeric\n"
                       ,"Unpredictable results will follow if casting to"
                       ,"numeric vector fails\n"))
    ppm <- as.numeric(ppm)
  }
  if (!is.numeric(Y)){
    cat(crayon::yellow("nmr.spectra.processing::calibrateSpectra >>"
                       ,"Non-numeric argument Y"))
    if (is.null(dim(Y))){
      cat(crayon::yellow(" begin cast as.numeric\n"
                         ,"Unpredictable results will follow if casting to"
                         ,"numeric vector fails"))
      Y <- as.numeric(Y)
    }
    else{
      cat(crayon::yellow(" being cast as.matrix\n"
                         ,"Unpredictable results will follow if casting to"
                         ,"numeric matrix fails\n"))
      Y <- as.matrix(Y)
    }
  }
  
  make.ref <- function(ref, frequency, j, cshift){
    if (ref == "tsp"){
      if (missing(cshift)) cshift <- 0
      return(new("NMRSignal1D"
                 ,id = "TSPbration"
                 ,chemicalShift = cshift
                 ,shape=list(name = "pseudoVoigt"
                             ,params=list(mu=.85,fwhm=0.0001,base=0))
                 ,peaks = list(
                   new("NMRPeak1D", x=cshift, y=1)
                 )
      ))
    }
    if (ref == "alanine"){
      if (missing(j)) j <- 7.26
      j <- j / (2 * frequency)
      if (missing(cshift)) cshift <- 1.48
      return(new("NMRSignal1D"
                 ,id = "Alabration"
                 ,chemicalShift = cshift
                 ,shape=list(name = "pseudoVoigt"
                             , params=list(mu=.85,fwhm=.6/frequency,base=0))
                 ,peaks = list(
                   new("NMRPeak1D", x=cshift - j, y=1)
                   ,new("NMRPeak1D", x=cshift + j, y=1)
                 )
      ))
    }
    if (ref %in% c("glucose","serum")){
     if (missing(j)) j <- 3.63
     j <- j / (2 * frequency)
     if(missing(cshift)) cshift <- 5.223
      return(new("NMRSignal1D"
                 ,id = "Glucobration"
                 ,chemicalShift = cshift
                 ,shape=list(name = "pseudoVoigt"
                             , params=list(mu=.85,fwhm=.6/frequency,base=0))
                 ,peaks = list(
                   new("NMRPeak1D", x=cshift - j, y=1)
                   ,new("NMRPeak1D", x=cshift + j, y=1)
                 )
      ))
    }
    cat(crayon::red("nmr.spectra.processing::calibrateSpectra >>"
                    ,"invalid reference\n"))
    stop()
  }
  
  if(missing(rOref)) rOref <- NULL
  
  #Assign tailored rOref if necessary
  if (is.null(rOref)){
    if (ref == "tsp"){
      rOref <- c(-0.02,0.02)
    }  else{
      if (ref == "serum"){
        rOref <- c(5.2,5.4)
      }
    }
  }
  
  #create ref if needed
  if (is.character(ref)){
    ref <- make.ref(ref, frequency, j, cshift)
  } 
  if (!("NMRSignal1D" %in% is(ref))){
    cat(crayon::red("nmr.spectra.processing::calibrateSpectra >>"
                    ,"invalid reference\n"))
    stop()
  }
  
  #Compute default rOref if necessary
  if(is.null(rOref)) rOref <- signalDomain(ref,30)
  
  #Align normalized spectra on rOref to reference and get shifts
  shifts <- calibrateToSignal(ppm,Y,ref,rOref, maxShift=maxShift
                              ,threshold=threshold,...)
  
  #Shift whole spectra by the corresponding shifts
  if (is.matrix(Y)){
    Y <- t(sapply(1:dim(Y)[1],function(i){
      shift <- shifts[i]
      y <- Y[i,]
      shiftSeries(y,shift, padding=padding, from=from)
    }))
    rownames(Y) <- rowLabels
    return(Y)
  } else{
    shiftSeries(Y,shifts,padding=padding,from=from)
  }
}

#Auxiliary function, returns the shift values but does not shift so that it can
#be used by the actual exports to either calibrate the spectra to a signal or
#calibrate the signal to a spectrum
calibrateToSignal <- function(ppm, Y, signal, rOref=signalDomain(signal,30)
                              , maxShift=1/3,threshold=0.2
                              ,lambda=1e3, ...){
  rOref <- crop(ppm,roi=rOref)
  ppm <- ppm[rOref]
  #bcorr annoyting broad features
  #Align normalized spectra on rOref to reference and get shifts
  if (is.matrix(Y)){
    Y <- Y[,rOref]
    Y <- baselineCorrection(Y,lambda=lambda)
    normS <- t(apply(Y,1,function(y) y / max(y)))
  } else{
    Y <- Y[rOref]
    Y <-  baselineCorrection(Y,lambda=lambda)
    normS <- Y / max(Y)
  }
  alignSeries(normS, signalToY(normalizeSignal(signal),ppm),shift=FALSE
              ,lag.max=length(ppm) * maxShift, threshold=threshold, ...)
}