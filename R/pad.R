pad <- function(x, n, side = 0, method="sampling"
                ,using=(length(x)*14/15):length(x),from=x){
  N <- length(x)
  if (!side %in% c(-1,0,1)){
    cat(crayon::red("nmr-spectra-processing::pad >>"
                    ,"wrong side specification: use -1 for left,"
                    ,"1 for right, and 0 for both"))
    stop()
  }
  if (method == "zeroes"){
    return(switch((side+2)
                  ,c(rep(0,n),x)
                  ,c(rep(0,n),x,rep(0,n))
                  ,c(x,rep(0,n))
    ))
  }
  if (method == "sampling"){
    #qc from
    if (!is.numeric(from)){
      cat(crayon::yellow("pad >>", "invalid argument value for 'from'",
                         "switching to default 'x' for sampling"))
      from <- x
    }
    #qc using
    if (is.logical(using)){
      if (length(using) != length(from)){
        cat(crayon::yellow("pad >>", "invalid argument value for 'using'",
                        "switching to default last 1/15 points for sampling"))
        using <- as.integer((length(from)*14/15):length(from))
      }
    }
    else{
      if (is.numeric(using)){
        using <- as.integer(using)
      }
      else{
        cat(crayon::yellow("pad >>", "invalid argument value for 'using'",
                           "switching to default c(9.5,10) for sampling"))
        using <- as.integer((length(from)*14/15):length(from))
      }
    }
    return(switch(side+2
                  ,c(sample(from[using],n,replace=TRUE)/2, x)
                  ,c(sample(from[using],n,replace=TRUE)/2, x 
                     , sample(from[using],n,replace=TRUE)/2)
                  ,c(x,sample(from[using],n,replace=TRUE)/2)
    ))
  }
  if (method == "circular"){
    #circular padding with more elements than the series contains makes no sense
    if (n>N){
      cat(crayon::red("nmr-spectra-processing::pad >>"
                      ,"n greater than series length not allowed"
                      ,"in circular padding"))
      stop()
    }
    #circular padding on both ends makes no sense
    if (side==-1) return(c(x[N-n+(1:n)],x))
    return(c(x,x[1:n]))
  }
}