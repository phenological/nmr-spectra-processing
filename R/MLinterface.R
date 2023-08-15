#Initialization of the global v8 context within the package namespace
#Necessary to interact with js and use ml-nmrSpectraProcessing capabilities
nsp.env <- new.env()
nsp.env$ct <- V8::v8()

#' Runs JS script on x-y data and returns the result on the R-side
#' 
#' xyApplyJS parses two R vectors into a ml xy object {x: [], y:[]}, passes the 
#' object to a JS environment using the V8 package, runs the given script with 
#' the xy object as input, captures the output and parses it back to R.
#'  Communication between R and JS is mediated by a V8 context, which means 
#'  that parsing is done through jsonlite.
#'  
#' @param x A vector (probably numeric), any other object that parses to a 
#' JS array through jsonlite::toJSON
#' @param y A vector (probably numeric) or any other object that parses to a 
#' JS array through jsonlite::toJSON
#' @param script A string of JS code. Within this string, "xy" represents the 
#' xy object built from the x and y parameters
#' @param v8Context The V8 context to be used to run the JS side. You should
#' not need to change this.
#' @returns The output of the script, parsed through jsonlite::fromJSON.
#' As a side effect, the variable xy will be (re)set on thev8Context
#' @examples xyApplyJS(x=c(1,2,3),y=c(1,2,3),script="{x:xy.x*2,y:xy.y-3}")
#' @import V8

xyApplyJS <- function(x,y,script){
  nsp.env$ct$assign("xy",list(x=x,y=y))
  nsp.evn$ct$get(JS(script))
}

mlPeakToPeak <- function(mlPeak){
  
}

mlSignalToSignal <- function(mlSignal){
  
}