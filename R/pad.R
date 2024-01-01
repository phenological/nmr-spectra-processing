#' Pads a series on either extreme with the number of specified "blank" points
#' 
#' @param x numeric vector, series to be padded
#' @param n numeric, number of points to be added
#' @param side numeric, whether to pad at the start (-1), at the end (1, default)
#'  or at both extremes (0)
#' @param method character, padding method. "zeroes" (default) pads with 0.
#'  "circular" pads the end with points from the start and the start with points
#'  from the end. "sampling" pads with points sampled from the input series
#' @param from logical or integer, optional. Filter selecting the region from x
#'  to be used in the "sampling" method. Default: the last 1/15th points
#' @returns the padded series
#' @export
pad <- function(x, n, side = c(-1,0,1)[3]
                ,method=c("zeroes","circular","sampling")[1]
                ,from=as.integer(length(x)*14/15):length(x)){
  if (!is.numeric(x)){
    cat(crayon::yellow("nmr-spectra-processing::pad >>"
                    ,"Argument x being cast as.numeric\n"
                    ,"Unpredictable results will follow if casting to"
                    ,"numeric vector fails\n"))
    x <- as.numeric(x)
  }
  if (!side %in% c(-1,0,1)){
    cat(crayon::red("nmr-spectra-processing::pad >>"
                    ,"wrong side specification: use -1 for left,"
                    ,"1 for right, and 0 for both\n"))
    stop()
  }
  n <- as.integer(n[1])
  if (method == "zeroes"){
    return(switch((side+2)
                  ,c(rep(0,n),x)
                  ,c(rep(0,n),x,rep(0,n))
                  ,c(x,rep(0,n))
    ))
  }
  if (method == "sampling"){
    # #qc from
    # if (!is.numeric(from)){
    #   cat(crayon::yellow("pad >>", "invalid argument value for 'from'",
    #                      "switching to default 'x' for sampling"))
    #   from <- x
    # }
    #qc using
    if (is.logical(from)){
      if (length(from) != length(x)){
        cat(crayon::yellow("pad >>", "invalid argument value for 'from'",
                        "switching to default last 1/15 points for sampling"))
        from <- as.integer(length(x)*14/15):length(x)
      }
    }
    else{
      if (is.numeric(from)){
        from <- as.integer(from)
      }
      else{
        cat(crayon::yellow("pad >>", "invalid argument value for 'using'",
                           "switching to default c(9.5,10) for sampling"))
        from <- as.integer(length(x)*14/15):length(x)
      }
    }
    return(switch(side+2
                  ,c(sample(x[from],n,replace=TRUE), x)
                  ,c(sample(x[from],n,replace=TRUE), x 
                     , sample(x[from],n,replace=TRUE))
                  ,c(x,sample(x[from],n,replace=TRUE))
    ))
  }
  if (method == "circular"){
    N <- length(x)
    #circular padding with more elements than the series contains makes no sense
    if (n>N){
      cat(crayon::red("nmr-spectra-processing::pad >>"
                      ,"n greater than series length not allowed"
                      ,"in circular padding"))
      stop()
    }
    #circular padding on both ends makes no sense; right padding by default
    if (side==-1) return(c(x[N-n+(1:n)],x))
    return(c(x,x[1:n]))
  }
}