#' Spectra calibration
#' Calibrates spectra
#' @param x numeric, spectra scale e.g. ppm
#' @param Y matrix, intensities, spectra in rows
#' @param ref character, reference signal for calibration. Supported references: 'glucose', 'alanine', 'tsp'
#' @param rOref numeric, optional. Limits of the Region of
#' Reference withing which the reference signal for calibration will be aligned.
#' Defaults: 5.15 - 5.3 for glucose, 1.4 - 1.56 for alanine, and -0.2 - 0.2 for tsp.
#' @param cshift numeric, chemical shift where the reference signal should be in the output
#' Defaults to the center of the rOref
#' @param j numeric. For doublet references, the coupling constant of the doublet
#' Defaults: .0065 for glucose and .0125 for alanine
#' @param threshold numeric, minimum cross-correlation requiered for alignment
#' see alignSeries
#' @param ... additional arguments for alignSeries (see details)
#' @details crops the Region of Reference of the spectra matrix and passes it
#' to alignSeries, along with a model of the reference signal, to calculate the
#' shifts that align the spectra on the rOref to the ref signal.
#' Then, it shifts the full spectra by the corresponding amounts.
#' @returns calibrated spectra matrix
#' @importFrom stats ccf
calibrateSpectra <- function(x, Y, ref, rOref, cshift, j, maxShift=1/3, threshold=0.2, ...){
  #standards for refs
  rOrefs <- list(glucose=c(5.15,5.3)
                 ,alanine=c(1.4,1.56)
                 ,tsp=c(-0.2,0.2)
                 )
  js <- list(glucose=0.0065,alanine=0.0125)
  
  #qc ref
  if (missing(ref)){
    cat(crayon::red("calibrateSpectra >>","ref must be"
                    ,"either 'glucose','alanine' or 'tsp'"))
    stop()
  }
  else{
    if (!(ref %in% c("alanine","glucose","tsp"))){
      cat(crayon::red("calibrateSpectra >>","ref must be"
                      ,"either 'glucose','alanine' or 'tsp'"))
      stop()
    }
  }
  
  #qc or get rOref
  if (missing(rOref)) rOref <- rOrefs[[ref]]
  else{
    if(!(is.numeric(rOref) & length(rOref==2))){
      cat(crayon::red("alignSeries >>",
                      "Invalid value for optional argument 'rOref':"
                      ,"must be a chem. shift. interval\n"
                      ))
        stop()
      }
    }

  #qc or compute cshift
  if (missing(cshift)) cshift <- mean(rOref)
  else{
    if (!is.numeric(cshift) | (cshift < min(x[rOref])) | (cshift > max(x[rOref]))){
      cat(crayon::red("calibrateSpectra >>", "Invalid cshift value"))
    }
  }
  
  #qc j
  if (missing(j)) j <- js[[ref]]
  else{
    if (!is.numeric(j)){
      cat(crayon::red("calibrateSpectra >>", "Invalid coupling constant j"))
    }
  }
  
  #make rOref filter
  rOref <- x >= rOref[1] & x <= rOref[2]
  
  #generate model reference peak
  #TBC: is it ok to always model with a lorentzian?
  if (ref=="tps") means <- cshift
  else means <- c(cshift - j/2, cshift + j/2)
  ref <- rowSums(lorentzian(x[rOref],mean=means
                   ,fwhm=0.001))
  
  #Align normalized spectra on rOref to reference and get shifts
  normS <- t(apply(Y[,rOref],1,function(y) y / max(y)))
  shifts <- alignSeries(normS, ref, shift=FALSE, lag.max=length(ref) * maxShift, threshold, ... )
  
  #Shift whole spectra by the corresponding shifts
  t(sapply(1:dim(Y)[1],function(i){
    shift <- shifts[i]
    y <- Y[i,]
    shiftSeries(y,shift)
  }))
}
