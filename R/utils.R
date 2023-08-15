

#' Create a filter to crop a spectral region
#' 
#' Quality of life function to spare you a few uncomfortable key strokes and some neural pulses. Questionable value.
#' @param x numeric, chemical shift scale e.g. ppm
#' @param roi numeric, optional. Upper and lower limit of the Region Of Interest (roi) to be cropped.
#' @param start, numeric, optional. Lower limit of the roi
#' @param end, numeric, optional. Upper limit of the roi
#' @details Either start, end or roi is required. Argument "roi" has priority. If start is given but not end,
#'  the lower limit is effectively set to min(x). If end is given but no start, the upper limit is effectively
#'  set to max(x).
#' @return logic, TRUE for elements of x within roi, FALSE otherwise
#' @export
crop <- function(x,start=-Inf,end=Inf,roi){
  if (missing(roi))
    return(x >= start & x <= end)
  return(x >= roi[1] & x <= roi[2])
}

#' Get the index of a chemical shift
#' 
#' Seeks the element of a discrete vector representing a continuous scale that is closest to a given value. Useful e.g. to get the approximate intensity of spectra at a given chemical shift.
#' Quality of life function to save you a few key strokes and some neural pulses. Questionable value.
#' @param x, numeric, chemical shift scale e.g. ppm
#' @param v, numeric, chemical shift value
#' @return integer, the index of the element of x that is closest to v
#' @export
getIdx <- function(x,v){
  which.min(abs(x-v))
}

#' Get the n spectra with the highest intensity on the given chemical shift range or value
#' 
#' ppm, numeric, chemical shift scale
#' Y, matrix, numeric, spectra intensities, spectra in rows, chemical shifts in columns
#' cshift, numeric, optional, chemical shift value to be maximized
#' n, integer, number of spectra to be returned, 10 by default
#' roi, numeric, optional, length 2 vector with the limits of the chemical shift range to be maximized
#' bottom, logic, if TRUE returns the n spectra with the lowest intensity
#' index, logic, whether to return the row indices of top spectra or the spectra themselves (see details)
#' @details If both a precise chemical shift (`cshift`) and a chemical shift range (`roi`) are passed, cshift takes priority
#' @return If index is TRUE, returns a vector with the row indices of the spectra with the highest (or lowest, see `bottom`) intensity at the given `cshift `or within the given `roi`. Otherwise (default) returns a matrix with the corresponding rows.
#' @export
top <- function(ppm,Y,cshift,n=10L,roi=c(-Inf,Inf),bottom=FALSE,index=FALSE){
  if (missing(cshift)){
    fi <- crop(ppm,roi=roi)
    idx <- order(apply(Y[,fi],1,max),decreasing = !bottom)[1:n]
  }
  else{
    idx <- order(Y[,getIdx(ppm,cshift)],decreasing = !bottom)[1:n]
  }
  if (index) idx else Y[idx,]
}

#' Fast spectra overlay plot
#' 
#' Wrapper to graphics::matplot for plotting spectra. Crude compared to ggplot2 but faster.
#' @param x numeric, your x scale e.g. ppm
#' @param y numeric or matrix, spectra intensities, one spectrum per row
#' @param roi, numeric, optional, limits of the Region Of Interest to be plotted. Defaults to the range of x.
#' @param by, numeric, optional, number of spectra to be overlayed on each plot
#' @param type, optional, defaults to "l" (lines), default recommended for speed
#' @param lty, optional, defaults to 1 (continuous line)
#' @param 
#' @param ..., additional arguments to be passed to matplot
#' @returns NULL
#' @import RColorBrewer
#' @export
smatplot <- function(x, y, by, roi, type="l",lty=1
                     ,reverse=TRUE,legend,label,palette="Set1",...){
  if (is.matrix(y))
    y <- t(y)
  else
    y <- as.matrix(y)
  
  #make color palette
  n <- RColorBrewer::brewer.pal.info[palette,"maxcolors"]
  colores <- RColorBrewer::brewer.pal(n=n,name=palette)
  
  if (!missing(roi)){
    fi <- x >= roi[1] & x <= roi[2]
    x <- x[fi]
    y <- y[fi,]
  }
  else{
    roi <- range(x)#c(min(x),max(x))
  }
  if (reverse){
    roi <- rev(roi)
  }
  
  if (missing(by)){
    matplot(x,y,type=type,lty=lty,xlim=roi,col=colores,...)
    if (!missing(legend)){
      if (missing(label))
        legend(legend,legend=1:dim(y)[2],text.col=colores)
      else
        legend(legend,legend=label[1:dim(y)[2]]
               ,text.col=colores)
    }
    
  } 
  else{
    n <- floor(dim(y)[2] / by)
    r <- dim(y)[2] %% by
    for (j in 1:n){
      soi <- 1:by + by*(j-1)
      matplot(x,y[,soi,drop=FALSE],type=type,lty=lty,col=colores,xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          legend(legend,legend=soi,text.col=colores)
        else
          legend(legend,legend=label[soi]
                 ,text.col=colores)
      }
    } 
    if (r>0){
      soi <- (by*n+1):(by*n+r)
      matplot(x,y[,soi,drop=FALSE],type=type,lty=lty,col=colores,xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          legend(legend,legend=soi,text.col=colores)
        else
          legend(legend,legend=label[soi],text.col=colores)
      }
    } 
  }
}

#' Calculate the domain of a NMRSignal1D
#' Estimates the domain of the signal as the coordinates of the outer peaks offset by n linewidths
#' @param aSignal, NMRSignal1D
#' @param n, numeric, linewidth multiplier used to offset the outer peaks
#' @returns numeric, the start and end of the signal's domain
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
#' @export
signalDomain <- function(aSignal, n=3){
  if ("NMRSignal1D" %in% is(aSignal)){
    offset <- aSignal@shape$params$fwhm * n
    offset <- c(-offset,offset)
    peaksRange <- range(sapply(aSignal@peaks, function(aPeak) aPeak@x))
    peaksRange + offset
  }
}

#WARNING: not exported, watch for collisions nevertheless (e.g. nmr-spectra-quantification::pseudoVoigt)
#' Gaussian function
#' 
#' Interpolates a gaussian with the given parameters on the given points.
#' Vectorized.
#' 
#' @export
#' @param x numeric, points to be interpolated
#' @param Max numeric, maxima of the gaussian(s)
#' @param Mean numeric, mean(s) of the gaussian(s)
#' @param fwhm numeric, full-width-at-half-max(s) of the gaussians(s)
#' @returns numeric matrix with gaussian(x) for each combination of parameters
#' @examples gaussian(1:20,Max=1,Mean=c(10,15))
gaussian <- function(x,max=1,mean=0,fwhm=1){
  mapply(function(max.,mean.,fwhm.) max. * exp(-4*log(2)*((x-mean.)^2)/fwhm.^2)
         , max, mean, fwhm)
}

#' Lorentzian function

#' Interpolates a lorentzian with the given parameters on the given points. 
#' Vectorized
#' @export
#' @param x numeric, points to be interpolated
#' @param Max numeric, maxima of the lorentzian(s)
#' @param Mean numeric, mean(s) of the lorentzian(s)
#' @param fwhm numeric, full-width-at-half-max(s) of the lorentzians(s)
#' @returns numeric matrix with lorentzian(x) for each combination of parameters
#' @examples lorentzian(1:20,Max=1,Mean=c(10,15))
lorentzian <- function(x,max=1,mean=0,fwhm=1){
  gamma2 <- fwhm/2
  gamma2 <- gamma2^2
  mapply(function(max.,mean.,gamma2.) max. * gamma2./((x-mean.)^2+gamma2.)
         , max, mean, gamma2)
}

#' Pseudo-Voigt function
#' 
#' Interpolates a Voigt function, using the pseudo-Voigt approximation, 
#' on the given points with the given parameters. 
#' The pseudo-Voigt approximation approximates a Voigt function as a linear 
#' combination of a lorentzian and a gaussian. 
#' Vectorized.
#' @export
#' @param x numeric points to be interpolated
#' @param mean numeric, mean of the distribution
#' @param fwhm numeric, full width at half max of the distribution
#' @param mu numeric, <= 1. The mu parameter of the pseudo-Voigt approximation.
#' The coefficientes of the linear combination are mu for the lorentzian and
#' 1 - mu for the gaussian
#' @returns numeric matrix with gaussian(x) for each combination of parameters
#' @examples pseudoVoigt(1:20,Max=1,Mean=c(10,15),mu=c(0,0.5))
pseudoVoigt <- function(x,max=1,mean=0,fwhm=1,mu=0){
  # if (mu==0)
  #   return(gaussian(x,max=max,mean=mean,fwhm=fwhm))
  # if (mu==1)
  #   return(lorentzian(x,max=max,mean=mean,fwhm=fwhm))
  mapply(function(max.,mean.,fwhm.,mu.)
    mu.*lorentzian(x,max.,mean.,fwhm.) + (1-mu.)*gaussian(x,max.,mean.,fwhm.)
    , max, mean, fwhm, mu)
}