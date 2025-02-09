#' Create a filter to crop a spectral region
#' 
#' Quality of life function to spare you a few uncomfortable key strokes and
#'  some neural pulses. Questionable value.
#' @param ppm numeric, chemical shift scale
#' @param roi numeric, optional. Upper and lower limit of the region of interest
#'  to be cropped.
#' @param start, numeric, optional. Lower limit of the region of interest.
#' @param end, numeric, optional. Upper limit of the region of interest.
#' @details Either \code{start}, \code{end} or \code{roi} is required. Argument 
#' \code{roi} has priority. If \code{start} is given but not \code{end}, the 
#' upper limit is effectively set to \code{max(ppm)}. If \code{end} is given but
#'  not \code{start}, the lower limit is effectively set to \code{min(ppm)}.
#' @returns logic, a filter for the elements of \code{ppm} within the \code{roi}.
#' @export
crop <- function(ppm,start=-Inf,end=Inf,roi){
  if (missing(roi))
    return(ppm >= start & ppm <= end)
  return(ppm >= roi[1] & ppm <= roi[2])
}

#' Get the index of a chemical shift
#' 
#' Quality of life function to spare you a few key strokes and some neural pulses.
#' @param ppm, numeric, chemical shift scale
#' @param ..., numeric, query chemical shift values
#' @details Seeks the element(s) of the chemical shift scale closest to the given value(s).
#' Useful e.g. to get the approximate intensity of spectra at a given c. shift.
#' @returns integer, the indices of the elements of ppm that are closest to the queried values
#' @export
getI <- function(ppm,...){
  sapply(c(...),function(v) which.min(abs(ppm-v)))
}

#' Get the \emph{n} spectra with the highest intensity on the given chemical shift
#' range or value
#' 
#' @param ppm, numeric, chemical shift scale.
#' @param Y, matrix, numeric, intensities, spectra in rows.
#' @param cshift, numeric, optional, query chemical shift value.
#' @param n, integer, number of spectra to be returned, 10 by default.
#' @param roi, numeric, optional, length 2 vector with the limits of the query
#' chemical shift range.
#' @param bottom, logic, if TRUE returns the \code{n} spectra with the \emph{lowest}
#' intensity instead.
#' @param index, logic, determnies whether to return the row indices of top
#'  spectra (TRUE) or the spectra themselves (FALSE, default).
#' @details If both a precise chemical shift and a chemical shift range are 
#' passed, \code{cshift} takes priority
#' @returns If index is TRUE, returns a vector with the row indices of the
#' \emph{n} spectra with the highest (lowest if \code{bottom}) intensity at
#'  the given \code{cshift} or  within the given \code{roi}. Otherwise (default)
#'  returns a matrix with the corresponding spectra.
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
#' This is similar to \code{\link[graphics]{matplot}} with \code{interactive==FALSE}
#' but uses a slightly different interface, allows to split spectra in multiple
#' plots, handles \code{...} parameters more accurately, implements a more
#' advanced \code{resolution} system, etc.
#' 
#' @param ppm numeric, spectra ppm scale
#' @param y numeric vector or matrix, nmr intensities, spectra in rows. If neither
#' vector nor matrix, argument is cast as.matrix. Type compliance is recommended.
#' @param roi numeric, optional, limits of the Region of Interest to be plotted.
#'  Defaults to the range of \code{ppm}.
#' @param by numeric, optional, number of spectra to be overlayed on each plot
#' @param lty optional, line type, defaults to 1 (continuous line). Note that this
#' parameter is passed both to \code{\link[graphics]{matplot}} and \code{\link[graphics]{legend}},
#' so you always get lines in your legend (see ... below)
#' @param legend, optional, position of the legend as specified in
#' \code{\link[graphics]{legend}}
#' @param label, character, optional, series labels for the legend
#' @param reverse logic, optional, if TRUE (default) the x scale of the plot 
#' increases from right to left as it's customary in NMR spectroscopy
#' @param resolution, character. If "full" all data are passed to 
#' \code{\link[graphics]{matplot}}. If "dev", data are binned to fit the graphic
#'  device's resolution before plotting. By default \code{smatplot} chooses
#'  according to data size. See Details.
#' @param reduce, function used to compute bin values. See Details
#' @param col vector of series colors. The default is \emph{Set1} copied from 
#' \code{\link[RColorBrewer]{RColorBrewer}}. 
#' @param ..., additional arguments for customization. These arguments are
#' passed either to \code{\link[graphics]{matplot}} for plot customization, or
#' to  \code{\link[graphics]{legend}} for legend customization. Beware that some
#' of these parameters are shared, such as \code{ltw} and the non-optional
#' parameters \code{lty} and \code{col}. This means, for instance, that it is 
#' not possible to get a legend without lines or to unmatch legend linewidths 
#' from the plot's.
#' @returns NULL
#' @details
#' tl;dr: For spectra with less than ~1.2 million points the default method
#' renders fast and exact plots. Otherwise:
#' * the default method adds a quick pre-processing step to make the plot fast 
#' to render, at the cost of inaccurate peak heights
#' * \code{resolution="dev", reduce=max} pre-processes slower, renders as fast 
#' and gives accurate peak heights, but introduces ppm inaccuracies visible when 
#' zooming in
#' * \code{resolution="full"} is exact but very slow to render.
#' 
#' Further details follow.
#' 
#' When working with high resolution spectra it is likely that the resolution
#' of the spectrum is higher than the pixel resolution of the graphic device,
#' which places an unnecessary burden on the renderer. By default \code{smatplot}
#' bins the spectra matrix to match the pixel resolution of the active graphic 
#' device (as reported by \code{\link[grDevices]{dev.size}}) if the input spectra
#' contain more than 1.2 million points total. The user can force full resolution
#' or binned resolution by setting \code{resolution} to "full" or "dev" 
#' respectively.
#' 
#' By default, bin values are computed by sampling each spectrum at regular 
#' intervals. This method is fast, which is good for the visualization of large 
#' datasets, but it produces inaccurate estimates of peak height. Alternatively, 
#' the user may pass a function to the \code{reduce} parameter that will be used 
#' to compute the bins. In this case, each spectrum is partitioned with 
#' \code{\link[base]{split}} and the \code{reduce} function is applied to compute 
#' each bin's intensity value and ppm. This adds a expensive pre-processing
#' step that makes the plot slower to compute but the result renders much 
#' faster than \code{resolution="full"}. Beware that spectra rendered in this manner
#' are chemical shift-accurate only up to the current graphic resolution.
#' This means that inaccuracies in the frequency coordinate may appear
#' if you enhance the resolution of the graphic device after plotting e.g. by
#' enlarging the plot window in RStudio.
#' From experience,\code{reduce=max} produces the most accurate picture of peak 
#' intensities, while \code{reduce=mean} or \code{median} provide a compromise 
#' between intensity and chemical shift accuracy.
#' 
#' @importFrom graphics matplot
#' @importFrom grDevices dev.size
#' @export
smatplot <- function(ppm, y, roi, by,lty=1,legend,label
                     ,reverse=TRUE
                     ,resolution=c(
                       "full","dev"
                       ,ifelse(length(as.matrix(y)) > 1.2e6, "dev", "full")
                       )[3]
                     ,reduce
                     ,col=c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00"
                                ,"#FFFF33","#A65628","#F781BF","#999999")
                     ,...){
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
    roi <- sort(roi)
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
  
  #Local wrappers to ease optional argument forwarding
  lmatplot <- function(...,fill,border,angle,density,bty,bg
                       ,box.lwd,box.lty,box.col,pt.bg,pt.cex,pt.lwd
                       ,xjust,yjust,x.intersp,y.intersp,adj,text.width
                       ,text.col,text.font,merge,trace
                       ,plot,ncol,horiz,title,insert,xpd,title.col,title.adj
                       ,title.cex,title.font,seg.len){
    matplot(...)
  }
  llegend <- function(...,adj,ann,cex.axis,cex.lab,cex.main,cex.sub
                      ,col.axis,col.lab,col.main,col.sub,crt,err,family,fg,font
                      ,font.axis,font.lab,front.main,font.sub,lab,las,lend,ljoin
                      ,lmitre,mai,mar,mex,mgp,mkh,page,smo,srt,tck,xaxp,xaxs
                      ,xaxt,xpd,yaxp,yaxs,yaxt
                      ,type,main,sub,xlab,ylab,asp
                      ,xlim,ylim,log,add,verbose,frame
                      ){
    legend(...)
  }
    
  
  if (missing(by)){
    lmatplot(x=ppm,y=y,type="l",xlim=roi,lty=lty,col=col,...)
    if (!missing(legend)){
      if (missing(label))
        llegend(x=legend,legend=1:dim(y)[2],col=col,...)
      else
        llegend(x=legend,legend=label[1:dim(y)[2]],lty=lty
               ,col=col,...)
    }
    
  } 
  else{
    n <- floor(dim(y)[2] / by)
    r <- dim(y)[2] %% by
    for (j in 1:n){
      soi <- 1:by + by*(j-1)
      lmatplot(x=ppm,y=y[,soi,drop=FALSE],type="l",lty=lty,col=col,xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          llegend(x=legend,legend=soi,col=col,lty=lty,...)
        else
          llegend(x=legend,legend=label[soi],col=col,lty=lty,...)
      }
    } 
    if (r>0){
      soi <- (by*n+1):(by*n+r)
      lmatplot(ppm,y[,soi,drop=FALSE],type="l",lty=lty,col=col,xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          llegend(x=legend,legend=soi,col=col,lty=lty,...)
        else
          llegend(x=legend,legend=label[soi],col=col,lty=lty,...)
      }
    } 
  }
}

#' Normalize a \code{\linkS4class{NMRSignal1D}}
#' 
#' Scales signal height to a maximum of 1
#' @param signal a \code{\linkS4class{NMRSignal1D}}
#' @returns scaled \code{\linkS4class{NMRSignal1D}}
#' @import methods
#' @importClassesFrom nmr.peaks NMRPeak1D
#' @importClassesFrom nmr.peaks NMRSignal1D
#' @export
normalizeSignal <- function(signal){
  summit <- max(sapply(signal@peaks,function(aPeak) aPeak@y))
  signal@peaks <- lapply(signal@peaks, function(aPeak){
    aPeak@y <- aPeak@y / summit
    return(aPeak)
  })
  return(signal)
}


#Now should be impoted from nmr.peaks
# gaussian <- function(x,mean=0,max=1,fwhm=1){
#   max * exp(-4*log(2)*(((x - mean) / fwhm)^2))
# }
# lorentzian <- function(x,mean=0,max=1,fwhm=1){
#   gamma2 <- fwhm/2
#   gamma2 <- gamma2^ 2
#   max * gamma2/((x - mean)^2 + gamma2)
# }
# pseudoVoigt <- function(x,mean=0,max=1,fwhm=1,mu=0){
#   (1 - mu) * gaussian(x, mean, max, fwhm) + mu * lorentzian(x, mean, max, fwhm)
# }