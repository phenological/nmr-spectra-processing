#' Shifts a series by the given value
#' 
#' @param x numeric, the series to be shifted
#' @param shift numeric, the amount to be shifted. Negative shifts to the left.
#' @param padding character, the method to be used to fill the empty extremes of
#'  the shifted series. See \code{\link{pad}} for details. Default: "zeroes".
#' @param from logical or integer, optional. Filter selecting the region from x
#'  to be used in the "sampling" method, see \code{\link{pad}} for details.
#'  Default: the last 1/15th points.
#' @returns numeric, the shifted series
#' @export
shiftSeries <- function(x,shift,padding="sampling"
                        ,from=as.integer(length(x)*14/15):length(x)
                        ){
  if (!is.numeric(x)){
    cat(crayon::yellow("nmr-spectra-processing::pad >>"
                       ,"Argument x being cast as.numeric\n"
                       ,"Unpredictable results will follow if casting to"
                       ,"numeric vector fails\n"))
    x <- as.numeric(x)
  }
  direction <- sign(shift)
  shift <- abs(shift)
  padded <- pad(x,shift,-direction,method=padding,from=from)
  if (direction==1) return(padded[1:length(x)])
  return(padded[(shift+1):length(padded)])
}

#' Shifts spectra by the given frequency
#' 
#' @param ppm numeric, chemical shift scale
#' @param Y numeric matrix or vector, NMR intensities, spectra in rows
#' @param shift numeric, the frequency to shift each spectrum by. Recycled as
#' necessary, see details
#' @param hertz logical. If TRUE, \code{shift} is interpreted as an
#' absolute frequency shift in hertz. If FALSE (default), it is interpreted as
#' a relative ppm chemical shift
#' @param SF numeric, field strength in megahertz for each spectrum. Used to
#'  convert \code{shift} to ppm units if hertz=TRUE. Recycled as necessary, see 
#'  details
#' @param interpol character, the interpolation method to use, see details.
#'  See  \code{\link[signal]{interp1}} for accepted values. The default "spline"
#'  should works best in most cases.
#' @details Works by adding shift to ppm and then interpolating Y at ppm assuming
#' Y aligned with the shifted ppm scale. Arguments \code{shift} and \code{SF} are
#' recycled to match the number of rows of Y. No warnings are raised, be responsible.
#' @returns The shifted spectra
#' @importFrom signal interp1
#' @export
shiftSpectra <- function(ppm,Y,shift,hertz=FALSE,SF=600
                         ,interpol="spline"){
  #Auxiliary function, shifts one spectrum
  shiftSpectrum <- function(ppm,y,shift,interpol){
    interp1(ppm+shift,y,ppm,method=interpol)
  }
  if(!is.matrix(Y)){
    if (is.vector(Y)){
      #conver to ppm if requested
      if (hertz) shift <- shift[1] / SF[1]
      #do the thing
      return(shiftSpectrum(ppm,Y,shift,interpol))
    } else{
      #Cast if not complying to type
      cat(crayon::yellow("nmr-spectra-processing::smatplot >>"
                         ,"Argument Y being cast as.matrix.\n"
                         ,"Unpredictable results may follow if casting to"
                         ,"numeric matrix fails\n"))
      Y <- as.matrix(Y)
    }
  }
  
  #repeat shift and SF as many times as necessary to match number of spectra
  ny <- dim(Y)[1]
  ns <- length(shift)
  nf <- length(SF)
  
  n <- ny %/% ns
  m <- ny %% ns
  shift <- rep(shift,n+sign(m))
  shift <- shift[1:ny]
  
  n <- ny %/% nf
  m <- ny%% nf
  SF <- rep(SF,n+sign(m))
  SF <- SF[1:ny]
  
  #convert to ppm if requested
  if(hertz){
    shift <- shift / SF
  }
  
  #do the thing
  t(
    sapply(1:(dim(Y))[1], function(i){
      shiftSpectrum(ppm,Y[i,],shift[i],interpol=interpol)
    })
  )
}

