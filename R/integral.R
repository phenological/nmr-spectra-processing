#' Integrate a region or fitted signal
#' 
#' @param x numeric or FALSE (default). Points used for numeric approximations. Separation between neighbor points (delta) is assumed to be constant. If length(x)==1, it is interpreted as a delta.
#' @param y a numeric or FALSE (default). Either signal or y must be given.
#' @param signal list. A ml-nmrProcessing signal parsed into R through a V8 context. Run about.signal() for details. Either signal or y must be given.
#' @param method character or function, the method of integration. Currently implemented:  "fwhm": full-width-at-half-maximum times maximum approximation on the fitted  signal  "rectangle": sum of the areas of rectangles centered at each y point, with  width equal to the difference between x points  "sum": sum of y
#' @returns numeric, the value of the integral
#' @export
#There are major inconsistencies with what I take to be a signal. This is to
#a large extent a reflection of the inconsistency of jsonlite's parsing
#and the kind of irreproducible objects it can generate
#I am workarounding it everywhere, a robust solution is needed
integral <- function(x,y,signal,method="fwhm",frequency=400){
  if (missing(signal)){
    if (missing(y)){
      cat(crayon::red("nmr-spectra-processing::integral>>"
                      , "Need either signal or y\n"))
      stop()
    }
    if (method=="rect"){
      if (missing(x)){
        cat(crayon::red("nmr-spectra-processing::integral>>"
                        , "Cannot use 'rect' method without 'x' argument\n"))
        stop()
      }
      if (length(x) == 1) return(sum(y)*x)
      return(sum(y)*(x[2]-x[1]))
    }
    if (method=="sum") return(sum(y))
    cat(crayon::red("nmr-spectra-processing::integral>>"
                    , "No valid method provided\n"))
    stop()
  }
  if(!is.null(signal$peaks)) signal <- signal$peaks
  if(!is.data.frame(signal) & is.list(signal)){
    signal <- jsonlite::fromJSON(jsonlite::toJSON(signal))
  }
  if (method == "fwhm"){
    #I assume all peaks are of the same shape because otherwise it's a madhouse!
    s <- signal$shape$kind[1]
    w <- signal$shape$fwhm / frequency
    h <- signal$y
    if (s == "pseudoVoigt"){
      mu <- signal$shape$mu
      return(sum(mu * w * h) * pi / 2 + sum((1-mu) * w * h) * 1.064467)
    }
    if (s == "lorentzian")
      return( pi * sum(w * h) / 2)
    if (s == "gaussian")
      return(1.064467 * sum(w * h))
    cat(crayon::red("integral>>","Invalid signal shape",s,"\n"))
  }
  if (missing(x)){
    cat(crayon::red("nmr-spectra-processing::integral>>"
                    , "Need x to interpolate signal to integrate using method"
                    ,method,"\n"))
    stop()
  }
  y <- signalsToY(x,list(signal),frequency)
  if (method == "rect") return(sum(y)*(x[2]-x[1]))
  if (method=="sum") return(sum(y))
  if (is.function(method)) return(method(x,y,signal))
  cat(crayon::red("nmr-spectra-processing::integral>>"
                  ,"No valid method provided\n"))
  stop()
}