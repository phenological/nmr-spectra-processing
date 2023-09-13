#' Integrate NMR signals
#' 
#' @param nmr a NMRSignalModel or a list with named elements x and y
#' @param method character or function, see details
#' @details Either conmputes the integrals of a fusion::NMRSignal1D, of
#' each signal in a fusion::NMRSignalModel (through sapply on signalsOutput),
#'  or of a NMR region specified as list(x,y). If method is a function, 
#'  this function is used to compute the integrals; in the case of a list
#'  the elements x and y are passed as arguments. If method is a string, it is 
#'  interpreted as a specifier for a standar method.Methods implemented in this 
#'  way are: "sum" (total intensity, not yet implemented for non-list arguments),
#'  "rect" (Riemann sum approximation, not yet implemente for non-lst arguments)
#'   and "fwhm" (full-width at half-maximum approximation, not valid for list)
#'  @return numeric, the integrals of the signals in the NMRSignalModel or the 
#' NMR region
#' @importClassesFrom fusion NMRPeak1D
#' @importClassesFrom fusion NMRSignal1D
#' @impotClassesFrom fusion NMRSignalModel
#' @export
setGeneric("integral", function(nmr,method){
  standardGeneric("integral")
})

#' @rdname integral-methods
#' @aliases integral,NMRSignal1D,character-method
setMethod("integral"
          ,signature(nmr = "NMRSignal1D", method="character")
          ,definition = function(nmr, method) {
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

#' @rdname integral-methods
#' @aliases integral,NMRPeak1D,character-method
setMethod("integral"
          ,signature(nmr = "NMRPeak1D", method="character")
          ,definition = function(nmr, method) {
            if (grepl("^f",method)){
              mu <- nmr@shape$param$mu
              w <- nmr@shape$param$fwhm
              h <- nmr@shape$param$y
              return(sum(mu * w * h) * pi / 2 + sum((1-mu) * w * h) * 1.064467)
            }
            cat(crayon::red("nmrSpectraProcessing::integral >>"
                            ,"invalid integration method for NMRPeak1D object"))
          })
          
#' @rdname integral-methods
#' @aliases integral,NMRSignalModel,character-method
setMethod("integral"
          ,signature(nmr = "NMRSignalModel", method = "character")
          ,definition = function(nmr,method){
            if (grepl("^s",method)){
              y <- sapply(nmr@signalsOutput, function(signal)
                signalToY(signal,nmr@ppm)
              )
              return(colSums(y))
            }
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

#' @rdname integral-methods
#' @aliases integral,NMRPeak1D,function-method
setMethod("integral"
          ,signature(nmr="NMRPeak1D", method = "function")
          ,definition = function(nmr,method){
            method(nmr)
          })

#' @rdname integral-methods
#' @aliases integral,NMRSignal1D,function-method
setMethod("integral"
          ,signature(nmr="NMRSignal1D", method = "function")
          ,definition = function(nmr,method){
            method(nmr)
          })

#' @rdname integral-methods
#' @aliases integral,NMRSignalModel,function-method
setMethod("integral"
          ,signature(nmr="NMRSignalModel", method = "function")
          ,definition = function(nmr,method){
            sapply(nmr@signalsOutput,method)
          })

#' @rdname integral-methods
#' @aliases integral,list,character-method
setMethod("integral"
          ,signature(nmr = "list", method = "character")
          ,definition = function(nmr,method){
            if (grepl("^s",method))
              return(sum(nmr$y))
            if (grepl("^r",method))
              return(sum(nmr$y) * (nmr$x[2] - nmr$x[1]))
            cat(crayon::red("nmrSpectraProcessing::integral >>"
                            ,"invalid method"))
          })

#' @rdname integral-methods
#' @aliases integral,list,function-method
setMethod("integral"
          ,signature(nmr = "list", method = "function")
          ,definition = function(nmr, method){
            method(nmr$x,nmr$y)
          })


            
#     if(is.null(signal$shape)){
#       s <- signal$kind[[1]]
#       w <- signal$fwhm / frequency
#     }
#     else{
#       s <- signal$shape$kind[[1]]
#       w <- signal$shape$fwhm / frequency
#     }
#     h <- signal$y
#     if (s == "pseudoVoigt"){
#       if (is.null(signal$shape)) mu <- signal$mu
#       else mu <- signal$shape$mu
#       return(sum(mu * w * h) * pi / 2 + sum((1-mu) * w * h) * 1.064467)
#     }
#     if (s == "lorentzian")
#       return( pi * sum(w * h) / 2)
#     if (s == "gaussian")
#       return(1.064467 * sum(w * h))
#     cat(crayon::red("integral>>","Invalid signal shape",s,"\n"))
#   }
#   if (missing(x)){
#     cat(crayon::red("nmr-spectra-processing::integral>>"
#                     , "Need x to interpolate signal to integrate using method"
#                     ,method,"\n"))
#     stop()
#   }
#   y <- signalsToY(x,list(signal),frequency)
#   if (method == "rect") return(sum(y)*(x[2]-x[1]))
#   if (method=="sum") return(sum(y))
#   if (is.function(method)) return(method(x,y,signal))
#   cat(crayon::red("nmr-spectra-processing::integral>>"
#                   ,"No valid method provided\n"))
#   stop()
#   
#   })
# }