#' Estimates the noise level of the spectrum
#' 
#' @param ppm numeric, chemical shift scale of the spectrum.
#' @param y numeric, spectrum intensities.
#' @param level numeric, quantile to estimate noise, specified as a probability
#' as in \code{\link[stats]{quantile}}. Default: 0.99.
#' @param rOref numeric, limits of the Region of Reference used to estimate 
#' noise. Default: 9.8 - 10 ppm.
#' @details Takes a 'blank' region of the spectrum as reference and
#'  estimates noise as the \code{level}-probability quantile of intensities
#'  in that region.
#' @returns numeric, estimation of the spectrum's level of noise.
#' @importFrom stats quantile
#' @export
noiseLevel <- function(ppm,y,level=0.99,rOref=c(9.8,10)){
  if (!is.numeric(ppm)){
    cat(crayon::yellow("nmr.spectra.processing::noiseLevel >>"
                       ,"Non-numeric argument ppm being cast as.numeric\n"
                       ,"Unpredictable results will follow if casting to"
                       ,"numeric vector fails\n"))
    ppm <- as.numeric(ppm)
  }
  if(!is.numeric(rOref) | length(rOref) != 2){
    cat(crayon::red("nmr.spectra.processing::noiseLevel >>"
                    ,"Invalid rOref\n"))
    stop()
  }
  fi <- ppm >= rOref[1] & ppm <= rOref[2]
  if (is.vector(y)){
    if (!is.numeric(y)){
      cat(crayon::yellow("nmr.spectra.processing::noiseLevel >>"
                         ,"Non-numeric argument y being cast as.numeric\n"
                         ,"Unpredictable results will follow if casting to"
                         ,"numeric vector fails\n"))
      y <- as.numeric(y)
    }
    return(quantile(y[fi],level,names=FALSE))
  }
  
  if (!is.matrix(y)){
    cat(crayon::yellow("nmr.spectra.processing::noiseLevel >>"
                       ,"Argument y being cast as.matrix\n"
                       ,"Unpredictable results may follow if casting to"
                       ,"numeric matrix fails\n"))
    y <- as.matrix(y)
  }  
  if(!is.numeric(y)){
    cat(crayon::red("nmr.spectra.processing::alingSeries >>"
                    ,"Expected y to be a numeric vector or matrix \n"))
    stop()    
    }  
  return(apply(y[,fi],1,function(r){
    quantile(r,level,names=FALSE)
    }))
}