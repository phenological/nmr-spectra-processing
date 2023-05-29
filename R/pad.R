pad <- function(x, n, side = 0, method="sampling",using=length(x)/15){
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
    return(switch(side+2
                  ,c(sample(x[1:using],n,replace=TRUE), x)
                  ,c(sample(x[1:using],n,replace=TRUE), x 
                     , sample(x[N-using+(1:using)],n,replace=TRUE))
                  ,c(x,sample(x[N-using+(1:using)],n,replace=TRUE))
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