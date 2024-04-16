#' Create a filter to crop a spectral region
#' 
#' Quality of life function to spare you a few uncomfortable key strokes and
#'  some neural pulses. Questionable value.
#' @param ppm numeric, chemical shift scale e.g. ppm
#' @param roi numeric, optional. Upper and lower limit of the Region of Interest
#'  to be cropped.
#' @param start, numeric, optional. Lower limit of the \code{roi}.
#' @param end, numeric, optional. Upper limit of the \code{roi}.
#' @details Either \code{start}, \code{end} or \code{roi} is required. Argument 
#' \code{roi} has priority. If \code{start} is given but not \code{end}, the 
#' lower limit is effectively set to \code{min(ppm)}. If \code{end} is given but
#'  not \code{start}, the upper limit is effectively set to \code{max(ppm)}.
#' @returns logic, a filter for the elements of \code{ppm} within the \code{roi}.
#' @export
crop <- function(ppm,start=-Inf,end=Inf,roi){
  if (missing(roi))
    return(ppm >= start & ppm <= end)
  return(ppm >= roi[1] & ppm <= roi[2])
}

#' Get the index of a chemical shift
#' 
#' Seeks the element of the chemical shift scale closest to the given value.
#' Useful e.g. to get the approximate intensity of spectra at a given c. shift.
#' Quality of life function to save you a few key strokes and some neural pulses.
#' Questionable value.
#' @param ppm, numeric, chemical shift scale
#' @param ..., numeric, chemical shift values to be found
#' @returns integer, the index of the element of ppm that is closest to v
#' @export
getI <- function(ppm,...){
  sapply(c(...),function(v) which.min(abs(ppm-v)))
}

#' Get the \emph{n} spectra with the highest intensity on the given chemical shift
#' range or value
#' 
#' @param ppm, numeric, spectra chemical shift scale.
#' @param Y, matrix, numeric, intensities, spectra in rows.
#' @param cshift, numeric, optional, chemical shift value to be maximized.
#' @param n, integer, number of spectra to be returned, 10 by default.
#' @param roi, numeric, optional, length 2 vector with the limits of the chemical
#'  shift range to be maximized.
#' @param bottom, logic, if TRUE returns the \code{n} spectra with the \emph{lowest}
#' intensity instead.
#' @param index, logic, whether to return the row indices of top spectra (TRUE)
#' or the spectra themselves (FALSE, default).
#' @details If both a precise chemical shift and a chemical shift range are 
#' passed, \code{cshift} takes priority
#' @returns If index is TRUE, returns a vector with the row indices of the spectra
#'  with the highest (or lowest, if \code{bottom}) intensity at the given 
#'  \code{cshift} or  within the given \code{roi}. Otherwise (default) returns a
#'   matrix with the corresponding spectra.
#' @export
top <- function(ppm,Y,cshift,n=10L,roi=c(-Inf,Inf),bottom=FALSE,index=FALSE){
  if (!is.numeric(ppm)){
    cat(crayon::yellow("nmr-spectra-processing::pad >>"
                       ,"Argument ppm being cast as.numeric\n"
                       ,"Unpredictable results may follow if casting to"
                       ,"numeric vector fails\n"))
    ppm <- as.numeric(ppm)
  }
  if (!is.matrix(Y)){
    cat(crayon::yellow("nmr-spectra-processing::pad >>"
                       ,"Argument Y being cast as.matrix\n"
                       ,"Unpredictable results may follow if casting to"
                       ,"numeric matrix fails\n"))
    Y <- as.matrix(Y)
  }
  if (!is.numeric(Y)){
      cat(crayon::red("nmr-spectra-processing::pad >>"
                         ,"Expected Y to be a numeric matrix\n"))
      stop()
  }
  n <- min(n,dim(Y)[1])
  if (missing(cshift)){
    fi <- crop(ppm,roi=roi)
    idx <- order(apply(Y[,fi],1,max),decreasing = !bottom)[1:n]
  }
  else{
    idx <- order(Y[,getI(ppm,cshift)],decreasing = !bottom)[1:n]
  }
  if (index) idx else Y[idx,]
}

#' Fast spectra overlay plot
#' 
#' Wrapper to \code{\link[graphics]{matplot}} for plotting spectra. Crude 
#' compared to \code{ggplot2} but faster.
#' 
#' @param ppm numeric, spectra ppm scale
#' @param y numeric vector or matrix, nmr intensities, spectra in rows. If neither
#' vector nor matrix, argument is cast as.matrix. Type compliance is recommended.
#' @param roi numeric, optional, limits of the Region of Interest to be plotted.
#'  Defaults to the range of \code{ppm}.
#' @param by numeric, optional, number of spectra to be overlayed on each plot
#' @param type optional, defaults to "l" (lines), default recommended for speed
#' @param lty optional, defaults to 1 (continuous line)
#' @param legend, optional, position of the legend as specified in
#' \code{\link[graphics]{legend}}
#' @param label, character, optional, series labels
#' @param reverse logic, optional, if TRUE (default) the x scale of the plot 
#' increases from right to left as it's customary in NMR spectroscopy
#' @param resolution, character. If "full" all data are passed to 
#' \code{\link[graphics]{matplot}}. If "dev", data are binned to fit the graphic
#'  device's resolution. By default\code{smatplot} chooses according to data 
#'  size. See Details.
#' @param reduce, function used to compute bin values. See Details
#' @param palette vector of colors, equivalent to matplot(col). The default is
#' Set1 copied from RColorBrewer
#' @param ..., additional arguments to be passed to \code{\link[graphics]{matplot}}
#' @returns NULL
#' @details
#' When working with high resolution spectra it is likely that the resolution
#' of the spectrum is higher than the pixel resolution of the graphic device,
#' which places unnecessary burden on the renderer. By default \code{smatplot}
#' bins the spectra matrix to match the pixel resolution of the active graphic 
#' device (as reported by \code{\link[grDevices]{dev.size}}) if the input spectra
#' contain more than 1.2 million points total. The user can force full resolution
#' or binned resolution by setting \code{resolution} to "full" or "dev" 
#' respectively.
#' 
#' By default, bin values are computed by sampling each spectrum at regular 
#' intervals. This method is fast, which is good for the visualization of large 
#' datasets, but it is inaccurate in the intensity coordinate. Alternatively, 
#' the user may pass a function to \code{reduce} that will be used to compute the
#' bins. In this case, each spectrum is partitioned with  \code{\link[base]{split}}
#' and the \code{reduce} function is applied to compute each bin's intensity value
#' and ppm. Beware that this may lead to inaccuracies in the frequency coordinate
#' that become apparent if you enhance the resolution of the graphic device while 
#' the plot is in view. From experience, \code{reduce=max} produces the most 
#' accurate picture of spectra peak intensities, while \code{reduce=mean} or
#' \code{median} provide a good compromise between accurate intensities and ppm. 
#' Either of the three should be equally accurate on the frequency scale
#' as long as the pixel resolution is not enhanced after plotting.
#' 
#' Keep in mind that the custom \code{reduce} plot may be slower to compute
#' but it will render faster.
#' @importFrom graphics matplot
#' @importFrom grDevices dev.size
#' @export
smatplot <- function(ppm, y, roi, by, type="l",lty=1,legend,label
                     ,reverse=TRUE
                     ,resolution=c(
                       "full","dev"
                       ,ifelse(length(as.matrix(y)) > 1.2e6, "dev", "full")
                       )[3]
                     ,reduce
                     ,palette=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"
                                ,"#FFFF33","#A65628","#F781BF","#999999"),...){
  #Cast if not complying to type
   if(!is.matrix(y) & !is.vector(y)){
    cat(crayon::yellow("nmr-spectra-processing::smatplot >>"
                       ,"Argument Y being cast as.matrix.\n"
                       ,"Unpredictable results may follow if casting to"
                       ,"numeric matrix fails\n"))
    y <- as.matrix(y)
   }
  #transpose for matplot compability
  if (is.matrix(y)) y <- t(y) else y <- as.matrix(y)
  #roi filter
  if (!missing(roi)){
    fi <- ppm >= roi[1] & ppm <= roi[2]
    ppm <- ppm[fi]
    y <- y[fi,,drop=FALSE]
  }
  else{
    roi <- range(ppm)
  }
  #Adjust resolution
  ##Choose resolution if not provided explicitly
  pointz <- dim(y)[1]
  smpls <- dim(y)[2]
  # if (missing(resolution)){
  #   if ((pointz * smpls) > 1.2e6){
  #     resolution <- "dev"
  #   } else{
  #     resolution <- "full"
  #   }
  # }
  
  #Adjust resolution if required
  if (resolution=="dev"){
    pixels <- grDevices::dev.size(units="px")[1]
    if (pointz > pixels){
      cat(crayon::yellow("nmr.spectra.processing::smatplot >>"
                         ,"Binning spectra to fit the graphic device's resolution\n"
                         ,"If you want to keep full resolution set argument resolution='full'\n"
                         )
          )
      pointsPerPixel <- pointz %/% pixels
      if (missing(reduce)){
        fi <- 1:pointz %% pointsPerPixel == 0
        ppm <- ppm[fi]
        y <- y[fi,,drop=FALSE]  
      } else{
        cat(crayon::yellow("nmr.spectra.processing::smatplot >>"
                           ,"Applying your custom reduce function\n"
                           )
            )
        bins <- as.factor(1:pointz %/% pointsPerPixel)
        ppm <- unname(sapply(split(ppm,bins),reduce))
        y <- apply(y,2,function(v){
          unname(
            sapply(
              split(v,1:pointz %/% pointsPerPixel)
              ,reduce
            )
          )
        })
      }
    }
  }
  #Reverse scale
  if (reverse){
    roi <- rev(roi)
  }
  
  if (missing(by)){
    matplot(ppm,y,type=type,lty=lty,xlim=roi,col=palette,...)
    if (!missing(legend)){
      if (missing(label))
        legend(legend,legend=1:dim(y)[2],text.col=palette)
      else
        legend(legend,legend=label[1:dim(y)[2]]
               ,text.col=palette)
    }
    
  } 
  else{
    n <- floor(dim(y)[2] / by)
    r <- dim(y)[2] %% by
    for (j in 1:n){
      soi <- 1:by + by*(j-1)
      matplot(ppm,y[,soi,drop=FALSE],type=type,lty=lty,col=palette,xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          legend(legend,legend=soi,text.col=palette)
        else
          legend(legend,legend=label[soi]
                 ,text.col=palette)
      }
    } 
    if (r>0){
      soi <- (by*n+1):(by*n+r)
      matplot(ppm,y[,soi,drop=FALSE],type=type,lty=lty,col=palette,xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          legend(legend,legend=soi,text.col=palette)
        else
          legend(legend,legend=label[soi],text.col=palette)
      }
    } 
  }
}

#' Calculate the domain of a \code{\linkS4class{NMRSignal1D}}
#' 
#' Estimates the domain of the signal as the coordinates of the outer peaks
#'  offset by \emph{n} linewidths (strictly: \emph{n} times full-width at half-max)
#' @param signal, \code{\linkS4class{NMRSignal1D}}
#' @param n numeric, linewidth multiplier used to offset the outer peaks
#' @returns numeric, the extremes of the signal's domain
#' @import methods
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
#' @export
signalDomain <- function(signal, n=3){
  if ("NMRSignal1D" %in% is(signal)){
    offset <- signal@shape$params$fwhm * n * c(-1,1)
    peaksRange <- range(sapply(signal@peaks, function(aPeak) aPeak@x))
    return(peaksRange + offset)
  }
  cat(crayon::red("signalDomain >>"
                  ,"Argument is not a S4 fusion::NMRSignal1D object"))
  stop()
}

#' Normalize a \code{\linkS4class{NMRSignal1D}}
#' 
#' Scales signal height to a maximum of 1
#' @param signal a \code{\linkS4class{NMRSignal1D}}
#' @returns scaled \code{\linkS4class{NMRSignal1D}}
#' @import methods
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
#' @export
normalizeSignal <- function(signal){
  summit <- max(sapply(signal@peaks,function(aPeak) aPeak@y))
  signal@peaks <- lapply(signal@peaks, function(aPeak){
    aPeak@y <- aPeak@y / summit
    return(aPeak)
  })
  return(signal)
}


## Now should be impoted from nmr-peak-fitting#' Gaussian function
#'
#' Interpolates a gaussian with the given parameters on the given points.
#' @param x numeric, points to be interpolated
#' @param max numeric, maxima of the gaussian(s)
#' @param mean numeric, mean(s) of the gaussian(s)
#' @param fwhm numeric, full-width-at-half-max(s) of the gaussians(s)
#' @returns numeric, vector of interpolated values
#' @examples gaussian(1:100,max=1,mean=0)
#' @export
gaussian <- function(x,mean=0,max=1,fwhm=1){
  max * exp(-4*log(2)*(((x - mean) / fwhm)^2))
}

#' Lorentzian function

#' Interpolates a lorentzian with the given parameters on the given points.
#' @param x numeric, points to be interpolated
#' @param max numeric, maxima of the lorentzian(s)
#' @param mean numeric, mean(s) of the lorentzian(s)
#' @param fwhm numeric, full-width-at-half-max(s) of the lorentzians(s)
#' @returns numeric, vector of interpolated values
#' @examples lorentzian(1:20,max=1,mean=0)
#' @export
lorentzian <- function(x,mean=0,max=1,fwhm=1){
  gamma2 <- fwhm/2
  gamma2 <- gamma2^ 2
  max * gamma2/((x - mean)^2 + gamma2)
}

#' Pseudo-Voigt function
#'
#' Interpolates a Voigt function, using the pseudo-Voigt approximation,
#' on the given points with the given parameters.
#'  The pseudo-Voigt approximation approximates a Voigt function as a linear
#' combination of a lorentzian and a gaussian.
#' @param x numeric points to be interpolated
#' @param max numeric, maxima of the lorentzian(s)
#' @param mean numeric, mean of the distribution
#' @param fwhm numeric, full width at half max of the distribution
#' @param mu numeric, <= 1. The mu parameter of the pseudo-Voigt approximation.
#' The coefficients of the linear combination are mu for the lorentzian and
#' 1 - mu for the gaussian
#' @returns numeric, vector of interpolated values
#' @examples pseudoVoigt(1:100,max=1,mean=0,mu=0.8)
#' @export
pseudoVoigt <- function(x,mean=0,max=1,fwhm=1,mu=0){
  (1 - mu) * gaussian(x, mean, max, fwhm) + mu * lorentzian(x, mean, max, fwhm)
}