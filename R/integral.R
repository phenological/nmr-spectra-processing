#' Integrate a region or fitted signal
#' 
#' @param x numeric or FALSE (default). Points used for numeric approximations. Separation between neighbor points (delta) is assumed to be constant. If length(x)==1, it is interpreted as a delta.
#' @param y a numeric or FALSE (default). Either signal or y must be given.
#' @param signal list. A ml-nmrProcessing signal parsed into R through a V8 context. Run about.signal() for details. Either signal or y must be given.
#' @param method character or function, the method of integration. Currently implemented:  "fwhm": full-width-at-half-maximum times maximum approximation on the fitted  signal  "rectangle": sum of the areas of rectangles centered at each y point, with  width equal to the difference between x points  "sum": sum of y
#' @returns numeric, the value of the integral
#' @export
#TBD: base R integrate as method, ml-side integration, other integration methods
integral <- function(x,y,signal,method="fwhm",frequency=400){
  if (missing(signal)){
    if (missing(y)){
      cat(crayon::red("integral>>", "Neither fitted signal nor y provided\n"))
      stop()
    }
    if (method=="rect"){
      if (missing(x)){
        cat(crayon::red("integral>>", "Cannot use 'rect' method without 'x' argument\n"))
        stop()
      } 
      if (length(x) == 1) return(sum(y)*x)
      return(sum(y)*(x[2]-x[1]))
    }
    if (method=="sum") return(sum(y))
    cat(crayon::red("integral>>", "No valid method provided\n"))
    stop()
  }
  #trick to make sure signal is a data.frame even if integral received a list
  if(!is.data.frame(signal) & is.list(signal)){
    signals <- jsonlite::fromJSON(jsonlite::toJSON(list(signal)))
  }
  if (method == "fwhm") return(sum(signal$y * signal$shape$fwhm / frequency))
  if (method == "rect"){
    if (missing(x)){
      cat(crayon::red("integral>>", "Cannot use 'rect' method without 'x' argument\n"))
      stop()
    }
    y <- signalsToY(x,signals,frequency)
    return(sum(y)*(x[2]-x[1]))
  }
  if (method=="sum") return(sum(y))
  if (is.function(method)) return(method(x,y,signal))
  cat(crayon::red("integral>>","No valid method provided\n"))
  stop()
}