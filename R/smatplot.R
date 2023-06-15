#' Fast spectra overlay plot
#' 
#' Wrapper to matplot for plotting spectra. Crude compared to ggplot2 but a lot faster.
#' @param x numeric, your x scale e.g. ppm
#' @param Y matrix, spectra matrix, one spectrum per row
#' @param limits, numeric, optional, limits of the region to be plotted. Defaults to full spectra.
#' @param by, numeric, optional, number of spectra to be overlay on each plot
#' @param type, optional, defaults to "l" (lines), default strongly recommended
#' @param ..., additional options to be passed to matplot
#' @returns NULL
#' @import RColorBrewer
#' @export
#Candidate to utilities.R
smatplot <- function(x, Y, by, roi, type="l"
                     ,reverse=FALSE,legend,label,palette="Paired",...){
  if (!missing(roi)){
    fi <- x >= roi[1] & x <= roi[2]
    x <- x[fi]
    Y <- Y[,fi]
  }
  else{
    roi <- TRUE
  }
  if (reverse){
    roi <- rev(roi)
  }
    
  if (missing(by)) matplot(x,t(Y),type=type,xlim=roi,...)
  else{
    colores <- brewer.pal(n=by,name=palette)
    n <- floor(dim(Y)[1] / by)
    r <- dim(Y)[1] %% by
    for (j in 1:n){
      soi <- 1:by + by*(j-1)
      matplot(x,t(Y[soi,]),type=type,col=colores[1:by],xlim=roi,...)
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
      matplot(x,t(Y[soi,]),type=type,col=colores[1:by],xlim=roi,...)
      if (!missing(legend)){
        if (missing(label))
          legend(legend,legend=soi,text.col=colores[1:by])
        else
          legend(legend,legend=label[soi],text.col=colores[1:by])
      }
    } 
  }
}