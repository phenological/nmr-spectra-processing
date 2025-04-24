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

#' #' Normalize a \code{\linkS4class{NMRSignal1D}}
#' #' 
#' #' Scales signal height to a maximum of 1
#' #' @param signal a \code{\linkS4class{NMRSignal1D}}
#' #' @returns scaled \code{\linkS4class{NMRSignal1D}}
#' #' @import methods
#' #' @importClassesFrom nmr.peaks NMRPeak1D
#' #' @importClassesFrom nmr.peaks NMRSignal1D
#' #' @export
#' normalizeSignal <- function(signal){
#'   summit <- max(sapply(signal@peaks,function(aPeak) aPeak@y))
#'   signal@peaks <- lapply(signal@peaks, function(aPeak){
#'     aPeak@y <- aPeak@y / summit
#'     return(aPeak)
#'   })
#'   return(signal)
#' }


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