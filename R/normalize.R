# 'Probabilistic Quotient Normalization
#' 
#' Computes the PQN estimation of the dilution factor for a spectrum or matrix of spectra
#' @param X numeric vector or matrix, input NMR
#' @param ref function or numeric vector, reference for the PQN. Can be provided
#' as either a function to be called on the spectra matrix to compute the reference
#' spectrum, or a vector giving the reference spectrum directly. Defaults to the
#' median spectrum.
#' @param of logical or integer, a filter to select the rows of the spectra matrix
#' to be used in the calculation of the reference.
#' @details See doi:10.1021/ac051632c
#' @returns numeric vector with the PQN dilution factors. Note that to apply the
#' correction each spectrum must be _divided_ by the corresponding dilution factor.
#' @export

pqn <- function(X, ref=median, of=TRUE){
  if (!is.numeric(X)){
    cat(crayon::red("nmrSpectraProcessing::pqn >>"
                    ,"Non-numeric spectrum"
                    ))
    stop()
  }
  if (is.function(ref)){
    if (!is.matrix(X)){
      cat(crayon::red("nmrSpectraProcessing::pqn >>"
                      ,"refuses to pqn a single spectrum with a reference"
                      ,"extracted from the same spectrum.\n"
                      ,"Please contact the developer if you think that would"
                      ,"make sense.\n"))
      stop()
    }
    ref <- apply(X[of,],2,ref)
  }
  
  if (is.vector(ref) & is.numeric(ref)){
    n <- ifelse(is.matrix(X),dim(X)[2],length(X))
    if (length(ref) == n){
      if (is.matrix(X)) return(apply(X,1,function(x) median(x / ref)))
      return(median(X / ref))
    }
    cat(crayon::red("nmrSpectraProcessing::pqn >>"
                    ,"Reference length", length(ref), "does not match",
                    "spectrum length",n,".\n"
                    ))
    stop()
  }
  cat(crayon::red("nmrSpectraProcessing::pqn >>"
                  ,"Invalid reference.\n"
                  ))
  stop()
}

#' Spectrum normalization
#' 
#' Convinence wrapper for spectra normalization
#' @param X numeric vector or matrix, input NMR
#' @param method, function to be used for normalization. Defaults to \code{\link{pqn}}
#' @param ..., additional parameters to method
#' @details if \code{X} is a matrix, \code{method} is applied to each row to obtain
#' the normalization coefficientes. Then, each row of \code{X} is divided by its
#' corresponding coefficient. If \code{X} is a vector, this functions is equivalent
#' to simply X / method(X). Since this is but a convenience wrapper, \code{method=max}
#' (max height normalization) and \code{method=sum} (total area normalization)
#' are parsed internaly so that they are applied to each row of the matrix
#'  rather than to the whole matrix in order to obtain the desired result.
#'  Total area normalization can also be obtained through \code{method=rowSums}
#' @returns numeric vector with the PQN dilution factors. Note that to apply the
#' correction each spectrum must be _divided_ by the corresponding dilution factor.
#' @export
normalize <- function(X,method=pqn,...){
  if (!is.numeric(X)){
    cat(crayon::red("nmrSpectraProcessing::pqn >>"
                    ,"Non-numeric spectrum"))
    stop()
  }
  if (is.matrix(X)){
    if (identical(method,max)){
      method <- function(m) apply(m,1,max)
    } 
    if (identical(method,sum)) method <- rowSums
  }
  X / method(X)
}
