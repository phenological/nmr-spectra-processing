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
#Candidate to utilities.R
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
    roi <- c(min(x),max(x))
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
      matplot(x,y[,soi],type=type,lty=lty,col=colores[1:by],xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          legend(legend,legend=soi,text.col=colores[1:by])
        else
          legend(legend,legend=label[soi]
                 ,text.col=colores)
      }
    } 
    if (r>0){
      soi <- (by*n+1):(by*n+r)
      matplot(x,y[,soi],type=type,lty=lty,col=colores[1:by],xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          legend(legend,legend=soi,text.col=colores[1:by])
        else
          legend(legend,legend=label[soi],text.col=colores[1:by])
      }
    } 
  }
}