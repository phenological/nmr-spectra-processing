#' Basic tool to annotate signals in a 1D NMR spectrum
#' 
#' Plots a spectrum with signals labeled
#' @param ppm numeric, chemical shift scale
#' @param y numeric, NMR intensities
#' @param annotations list, text annotations. Can be expressions.
#' @param peaks list, coordinates of the peaks associated with each annotation
#' @param roi numeric, optional. Upper and lower limit of the Region of Interest
#' to be plotted
#' @param tcol, text color for annotation
#' @param lcol, line color for annotation
#' @param window, numeric, ppm window for automatic annotation position adjustement.
#'  Default 0.002. See details.
#' @param ylim, by default ("auto") plotting area is adjusted to make sure all
#' labels fit. You may override this behavior by providing explicit limits
#' as in \code{\link[base]{plot}}
#' @param add, logic, should the annotation be added to the active plot (TRUE) or a
#' new plot created (FALSE)? This is useful to manually position each individual label
#' for optimal visualization, to pass different graphic parameters (e.g. lwd) for
#' spectrum trace and for annotation, or to keep a cleaner code in extensive annotations. 
#' Default FALSE.
#' @param adj, numeric, positioning of labels, see \code{\link[graphics]{par}}
#' @param delta, horizontal shift, in chemical shift units, of the text labels relative
#'  to the peak positions
#' @param epsilon, numeric, scaling factor for marker length, scales it
#' proportionally to the plot window. More precisetly, this controls the vertical
#' spacing between the top peak and the text label, so the actual length
#' of the marker depends on other factors such as `delta` positioning. Default 0.02 
#' @param marker, character, style of the annotation marker. Default: "trinket"
#' @param ... additional graphical parameters, see \code{\link[graphics]{par}}, 
#' see \code{\link{smatplot}}.
#' @details Annotations are re positioned to the nearest maximum intensity within
#' \code{window} of the given chemical shifts. You may compare the default
#' with  \code{window=0} to see the effect. Since annotations need the spectrum
#'  for adjustment, arguments \code{ppm} and \code{y} are mandatory even if
#' \code{add=TRUE}. For a similar reason, you should pass the original
#'  \code{ylim} to further calls using \code{add=TRUE}.
#' @returns NULL
#' 
#' @importFrom graphics lines text
#' @export
sannotate <- function(ppm,y,annotations,peaks,roi,tcol="black",lcol="black",window=0.002
                      ,ylim="auto",add=FALSE,adj=c(0.5,-0.5)
                      ,delta=0,epsilon=0.02,marker=c("trinket","range")[1],...){
  #qc
  m <- length(annotations)
  if (m != length(peaks)){
    cat(crayon::red("nmr.spectra.processing::sannotate>>"
                    ,"length of peaks and annotations lists must match\\n"))
  }
  
  #Local wrappers to ease optional argument forwarding
  llines <- function(...,lcol,col,roi, by, lty, legend, label, reverse, pos
                     ,resolution, reduce,adj, offset, vfont,font
                     ,xlim,ylim){
    lines(...,col=lcol)
  }
  ltext <- function(...,tcol,col){
    text(...,col=tcol)
  }
  lsmatplot <- function(...,add,adj,offset,vfont,font, pos){
    smatplot(...)
  }
  
  # #Compute window in index scale
  # window <- as.integer(window / (ppm[2] - ppm[1]))
  
  #Cut to roi
  if (missing(roi)) roi <- range(ppm)
  fi <- crop(ppm,roi=roi)
  ppm <- ppm[fi]
  y <- y[fi]
  
  # m <- length(annotations)
  #Adjust peak coordinates to nearest maxima within the window provided
  for (i in 1:m){
    peaks[[i]] <- sapply(peaks[[i]], function(px){
      # xi <- getI(ppm,px)
      xmin <- getI(ppm,px-window)
      xmax <- getI(ppm,px+window)
      xi <- xmin:xmax#(xi - window):(xi + window)
      return(ppm[xi[which.max(y[xi])]])
    })
  }
  #Signal center (named "means" for historical reasons)
  sigXmeans <- sapply(peaks,function(x) mean(range(x)))
  #Signal maxima
  if(marker == "range"){
    sigYmaxs <- sapply(peaks, function(x) max(y[crop(ppm,roi=range(x))]))
  } else{
    sigYmaxs <-sapply(peaks, function(x) max(y[getI(ppm,x)]))
  }
  
  # #Signal labels
  # l <- names(ann)
  
  #Estimate ideal ymax to ensure annotation fits (unless ylim provided)
  if (ylim[1] == "auto"){
    ysupp <- max(max(y), max(sigYmaxs) / (1-5*epsilon)#max(sigYmaxs) + minLength * 5 + poffset
                 )
    ylim <- c(min(y),ysupp)
  }
  
  #Viz spacing parameters
  ref <- diff(ylim) #reference: vertical plot
  poffset <- ref * 0.01 #epsilon/2 #peak-to-marker separation
  minLength <- ref * epsilon #vertical marker length (half of it, actually)

  if (!add) lsmatplot(ppm,y,ylim=ylim,...)
  # #Initialize empty plot on the desired range
  # plot(range(ppm),c(min(y),max(y) + minLength * 5 + poffset),type="n"
  #      ,xlab="ppm",ylab="")
  
  #Add labels to the annotation
  ltext(x=sigXmeans + delta,y=sigYmaxs + 2*minLength + poffset,labels=annotations,tcol=tcol
       ,adj=adj,...)
  #Add markers
  for (i in 1:m){
    ta <- peaks[[i]]
    #"range" style
    if (marker == "range"){
      #Add horizontals
      llines(range(ta), rep(sigYmaxs[[i]] + minLength + poffset,2),lcol=lcol,...)
      #Add verticals
      llines(rep(min(ta),2),c(sigYmaxs[[i]] + poffset #+ minLength/2
                              ,sigYmaxs[[i]] + poffset + minLength*2#*3/2
                              )
             ,lcol=lcol,...
             )
      llines(rep(max(ta),2),c(sigYmaxs[[i]] + poffset #+ minLength/2
                              ,sigYmaxs[[i]] + poffset + minLength*2#3/2
      )
      ,lcol=lcol,...
      )
    } else{
      # "trinket" style
      #Trinkets for mult > 1
      if (length(ta) > 1){
        ##Add horizontals
        llines(range(ta), rep(sigYmaxs[[i]] + minLength + poffset,2),lcol=lcol,...)
        ##Add branches
        for (x in ta){
          llines(rep(x,2),c(y[getI(ppm,x)] + poffset
                            ,sigYmaxs[[i]] + minLength + poffset
          )
          ,lcol=lcol,...
          )
        }
        ##Add trunk
        # llines(rep(mean(ta),2),c(sigYmaxs[[i]] + minLength + poffset
        llines(c(sigXmeans[[i]],sigXmeans[[i]]+delta),c(sigYmaxs[[i]] + minLength + poffset
                                                        ,sigYmaxs[[i]] + 2 * minLength + poffset
        )
        , lcol=lcol,...
        )
      } else{
        #For singlets, a single vertical
        llines(c(ta,ta+delta),c(sigYmaxs[[i]] + poffset
                                ,sigYmaxs[[i]] + 2 * minLength + poffset
        )
        ,lcol = lcol,...
        )
      }
    }
  }
}

# eg <- list("Valine?"=c(1.013,1.0305,1.064,1.082)
#            ,"Alanine"=c(1.485,1.503)
#            ,"Ethanol"=c(1.1822,1.1996,1.217)
#            ,"Morcilline"=c(1.3055)
#            ,"Climberino"=c(1.336,1.351)
#            )
# 
# #Testing optimization of label position to avoid overlap with the spectrum
# ##Get label bounding box (in the current plot)
# bbox <- function(label, x, y, resolution){
#   h <- strheight(label)
#   w <- strwidth(label)
#   
# }

##Compute a descriptor of overlap based on cross-correlation of the NMR to the
##label's bounding box


