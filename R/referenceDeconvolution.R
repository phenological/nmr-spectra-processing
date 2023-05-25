#' @importFrom gsignal hilbert
#' @export
referenceDeconvolution <- function(x,y,rOref,frequency=400
                                   ,signals=list(list()),zero=.001,...){
  #inverse fourier transform with normalization
  ifft <- function(x){
  fft(x, inverse = TRUE) / length(x)
}
  #qc input lengths
  if(length(x)!=length(y)){
    cat(crayon::red("referenceDeconvolution >>", "lengths of x and y don't match\n"))
    stop()
  }
 
  #avoid zeros in y
  #y[y==0] <- zero
   
  #parse and construct rOref
  if (missing(rOref)) rOref = TRUE
  else{
    if (!is.logical(rOref)){
      if (is.numeric(rOref) & length(rOref==2))
        rOref <- x >= rOref[1] & x <= rOref[2]
      else{
        cat(crayon::red("referenceDeconvolution >>", "invalid rOref\n"))
        stop()
      }
    }
  }
  
  #calculate trailing zeros outside of the rOref
  trail_l <- rep(0,which.max(rOref) - 1)
  trail_r <- rep(0,which.max(rev(rOref)) - 1)
  
  #get reference, fill trailing zeroes
  ref <- c(trail_l, y[rOref], trail_r)
  
  #temporary workaround: if there are ppm<0 in the rOref, shift the spectrum 
  #to start at zero, fit the signal and then shift back
  hasNegative <- sum(x<0) > 0
  if (hasNegative){
    shift <- x[1]
    x <- x - shift
  }
  #fit model
  refModel <- fitSignals(x,y,signals=signals,roi=rOref,frequency=frequency,...)
  if (is.null(refModel)){
    cat(crayon::yellow("referenceDeconvolution>>"
                       ,"signal fitting with fitSignals in the rOref returned NULL\n"
                       ,"deconvolution aborted\n"
                       ,"returning the input spectrum\n"))
    return(y)
  }
  refModel <- signalsToY(x,refModel$signals,frequency=frequency)
  #shift the scale back to its original position, if necessary
  if (hasNegative) x <- x + shift
  
  #move to time domain
  ref <- ifft(hilbert(ref))
  refModel <- ifft(hilbert(refModel))
  y <- ifft(hilbert(y))
  
  #deconv., correct and go back to freq. domain
  Re(fft(y * refModel / ref))
}
