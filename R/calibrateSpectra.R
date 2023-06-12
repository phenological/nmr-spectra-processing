#' Spectra calibration
#' Calibrates spectra
#' @param x numeric, spectra scale e.g. ppm
#' @param Y matrix, intensities, spectra in rows
#' @param ref character, semi-optional. Reference signal for calibration. Supported references: 'glucose', 'alanine', 'tsp'. One of ref and cshift must be provided.
#' @param rOref numeric, optional. Limits of the Region of Reference within which the reference signal for calibration will be aligned. Defaults: 5.15 - 5.3 for glucose, 1.4 - 1.56 for alanine, and -0.2 - 0.2 for tsp.
#' @param cshift numeric, chemical shift where the reference signal should be in the output. Defaults to the center of the rOref
#' @param j numeric. For doublet references, the coupling constant of the doublet. Defaults: .0065 for glucose and .0125 for alanine. 
#' @param threshold numeric, minimum cross-correlation requiered for alignment. See alignSeries
#' @param method character, the method to be used to fill the empty extremes of the shifted series (padding). See 'pad' for details. Default: "sampling".
#' @param using numeric, the number of points from the extremes to be used for padding in the "sampling" method. See 'pad' for details. Default: 1/15th of the series' length.
#' @param plot logical. Argument passed to ccf through alignSeries; if TRUE, each spectrum's correlation to the reference is plotted. Default: FALSE. You should not need to change the default unless you are getting misalignments and you want to check the cross-correlation to fine-tune the alignment parameters.
#' @param ... additional arguments for alignSeries, see Details
#' @details crops the Region of Reference of the spectra matrix and passes it to alignSeries, along with a model of the reference signal, to calculate the shifts that align the spectra on the rOref to the ref signal. Then, it shifts the full spectra by the corresponding amounts.
#' @returns calibrated spectra matrix
#' @importFrom stats ccf
calibrateSpectra <- function(x, Y, ref, rOref, cshift, j, fwhm
                             , maxShift=1/3, threshold=0.2, padding="sampling"
                             , using=c(9.5,10), from=apply(Y,2,median)
                             , plot=FALSE,...){
  #parse using as chem. shift. Any other form will be left for pad to parse.
  if (is.numeric(using) & length(using) == 2)
    using <- crop(x,using)
  
  #standards for refs
  rOrefs <- list(glucose=c(5.15,5.3)
                 ,alanine=c(1.4,1.56)
                 ,tsp=c(-0.1,0.1)
                 )
  js <- list(glucose=0.0065,alanine=0.0125)
  
  #qc ref
  if (missing(ref)){
    if (missing(rOref)){
      cat(crayon::red("calibrateSpectra >>","you must provide a reference signal"
                      ,"for calibration; either a standard"
                      ,"ref=c('glucose','alanine', 'tsp'), or a custom."
                      ,"reference via rOref, cshift and j\n"))
      stop()
    }
  }
  else{
    if (!(ref %in% c("alanine","glucose","tsp"))){
    cat(crayon::red("calibrateSpectra >>","supported ref are:"
                    ,"'glucose','alanine' or 'tsp'\n"))
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
    if (!is.numeric(cshift) | (cshift < rOref[1]) | (cshift > rOref[2])){
      cat(crayon::red("calibrateSpectra >>", "Invalid cshift value\n"))
      print(is(cshift))
      stop()
    }
  }
  
  #qc j
  if (missing(j)){
    if (missing(ref)){
      j=NULL
    }
    else{
      if (ref %in% c("alanine","glucose")) j <- js[[ref]]
      else j=NULL
    }
  }
  else{
    if (!is.numeric(j) | length(j)!=1){
      cat(crayon::red("calibrateSpectra >>", "Invalid coupling constant j\n"))
      stop()
    }
  }
  
  #qc fwhm
  if (missing(fwhm)) fwhm=0.0017
  else{
    if (!is.numeric(fwhm) | length(fwhm)!=1){
      cat(crayon::red("calibrateSpectra >>", "Invalid fwhm\n"))
      stop()
    }
  }
  
  #make rOref filter
  rOref <- x >= rOref[1] & x <= rOref[2]
  
  #generate model reference peak
  #If no j, model a singulet
  #Warning: this is using a lorentzian? Good for now, but we should gravitate to
  #pseudoVoigt once the problems with the optimization have been solved
  if (is.null(j)) means <- cshift
  else means <- c(cshift - j/2, cshift + j/2)
  ref <- rowSums(lorentzian(x[rOref],mean=means
                   ,fwhm=fwhm))
  
  #Align normalized spectra on rOref to reference and get shifts
  normS <- t(apply(Y[,rOref],1,function(y) y / max(y)))
  shifts <- alignSeries(normS, ref, shift=FALSE
                        , lag.max=length(ref) * maxShift, threshold=threshold
                        , plot=plot, ...)
  
  #Shift whole spectra by the corresponding shifts
  t(sapply(1:dim(Y)[1],function(i){
    shift <- shifts[i]
    y <- Y[i,]
    shiftSeries(y,shift, padding=padding, using=using, from=from)
  }))
}
