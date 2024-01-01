#' Shifts a series by the given value
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