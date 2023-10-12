#' Integrate NMR signals
#' 
#' @param nmr a \code{\linkS4class{NMRSignalModel}}, \code{\linkS4class{NMRSignal1D}}
#'  or a list with named elements \emph{x} and \emph{y}
#' @param method character or function, integration method to be used. 
#' If a string, it is interpreted as a specifier for an implemented
#' method: 
#' \itemize{
#' \item "rect", Riemann sum approximation; 
#' \item "fwhm", full-width at half-maximum approximation (not valid for list)
#' }
#' If a function, this function is used to compute the integrals with \code{nmr}
#' as argument (in the case of a list the elements, x and y are passed as separate
#' arguments.
#' @param ... additional arguments
#' @details In the "rect" method, signals are constrained to either the domain
#' of the \code{\linkS4class{NMRSignalModel}}, or to 3 fwhm around the
#' \code{\linkS4class{NMRSignal1D}}'s outer peaks
#' 
#'  For \code{\linkS4class{NMRSignalModel}} and list(x,y), the ppm resolution
#'  is read from the argument. For \code{\linkS4class{NMRSignal1D}}, 
#'  a resolution of 0.00023 ppm is used.
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
#' @importClassesFrom fusion NMRSignalModel
# importFrom nmr.peak.fitting gaussian lorentzian pseudoVoigt
#' @export
setGeneric("integral", function(nmr,method,...){
  standardGeneric("integral")
})

#' @describeIn integral NMRSignal1D method
#' @param delta numeric, the chemical shift gap between interpolated points. 
#' Though the default should work in most scenarios, you may want to adjust it to
#' the resolution of your experimental spectra.
#' @param offset numeric, a factor of the signal's full-width at half-maximum
#' (fwhm) used to determine the domain of integration: the signal is integrated
#' in a domain that offsets its outer peaks by \code{offset * fwhm}.
#' @returns numeric, the signal's integral
setMethod("integral"
          ,signature(nmr = "NMRSignal1D", method="character")
          ,definition = function(nmr, method, delta=0.00023, offset=3) {
            if (grepl("r",method)){
              ppm <- signalDomain(nmr,offset)
              ppm <- seq(ppm[1],ppm[2],by=delta)
              return(sum(signalToY(nmr,ppm)) * delta)
            }
            if (grepl("^f",method)){
              mu <- nmr@shape$param$mu
              w <- nmr@shape$param$fwhm
              h <- sapply(nmr@peaks, function(p) p@y)
              return(sum(mu * w * h) * pi / 2 + sum((1-mu) * w * h) * 1.064467)
            }
            cat(crayon::red("nmrSpectraProcessing::integral >>"
                            ,"invalid integration method"
                            ,"for NMRSignal1D object"))
          })

# #' @rdname integral-methods
# #' @aliases integral,NMRPeak1D,character-method
# setMethod("integral"
#           ,signature(nmr = "NMRPeak1D", method="character")
#           ,definition = function(nmr, method) {
#             if (grepl("^f",method)){
#               mu <- nmr@shape$param$mu
#               w <- nmr@shape$param$fwhm
#               h <- nmr@shape$param$y
#               return(sum(mu * w * h) * pi / 2 + sum((1-mu) * w * h) * 1.064467)
#             }
#             cat(crayon::red("nmrSpectraProcessing::integral >>"
#                             ,"invalid integration method for NMRPeak1D object"))
#           })
          
#' @describeIn integral NMRSignalModel method
#' @returns numeric, a vector of integrals for each signal in the NMRSignalModel
setMethod("integral"
          ,signature(nmr = "NMRSignalModel", method = "character")
          ,definition = function(nmr,method){
            if (grepl("^r",method)){
              y <- sapply(nmr@signalsOutput, function(signal)
                signalToY(signal,nmr@ppm)
              )
              return(colSums(y) * (nmr@ppm[2] - nmr@ppm[1]))
            }
            if (grepl("^f",method)){
              return(
                sapply(nmr@signalsOutput, function(signal){
                  mu <- signal@shape$param$mu
                  w <- signal@shape$param$fwhm
                  h <- sapply(signal@peaks, function(p) p@y)
                  sum(mu * w * h) * pi / 2 + sum((1-mu) * w * h) * 1.064467
                })
              )
            }
            cat(crayon::red("nmrSpectraProcessing::integral >>"
                            ,"invalid integration method for"
                            ,"NMRSignalModel1D object"))
          })

# #' @rdname integral-methods
# #' @aliases integral,NMRPeak1D,function-method
# setMethod("integral"
#           ,signature(nmr="NMRPeak1D", method = "function")
#           ,definition = function(nmr,method){
#             method(nmr)
#           })

#' @describeIn integral NMRSignal1D method with custom integration function
setMethod("integral"
          ,signature(nmr="NMRSignal1D", method = "function")
          ,definition = function(nmr,method){
            method(nmr)
          })

#' @describeIn integral NMRSignalModel method with custom integration function
setMethod("integral"
          ,signature(nmr="NMRSignalModel", method = "function")
          ,definition = function(nmr,method){
            sapply(nmr@signalsOutput,method)
          })

#' @describeIn integral x,y data series method
#' @returns numeric, the integral of the data series
setMethod("integral"
          ,signature(nmr = "list", method = "character")
          ,definition = function(nmr,method){
            if (grepl("^r",method))
              return(sum(nmr$y) * (nmr$x[2] - nmr$x[1]))
            cat(crayon::red("nmrSpectraProcessing::integral >>"
                            ,"invalid method"))
          })

#' @describeIn integral x,y data series method with custom integration function
setMethod("integral"
          ,signature(nmr = "list", method = "function")
          ,definition = function(nmr, method){
            method(nmr$x,nmr$y)
          })