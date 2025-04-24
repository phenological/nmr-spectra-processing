#' Visualize spectra in nmrium
#' 
#' Quality of life wrapper to \code{hastsaLaVista2} utils
#' @param ppm numeric, chemical shift scale
#' @param yr, matrix, numeric, nmr intensities, real part, samples in rows
#' @param yi, matrix, numeric, nmr intensities, imaginary part, samples in rows
#' Default: NULL
#' @param roi numeric, optional. Upper and lower limit of the region of interest
#'  to be cropped. Default \code{range(ppm)}
#' @param cache, directory to cache the data for visualization. A http server 
#' will be setup in this directory. Default \code{getwd()}
#' @param metadata, data.frame, optional, additional experimental data to be 
#' exported. Samples in rows.
#' @param ..., additional arguments for server initialization and data export
#' (see \code{hastaLaVista2::exportReIm})
#' @details R is a powerful software for data analysis, however it is not so
#' good at interactive plotting (we don't like plotly). This function interfaces
#' with nmrium to allow interactive spectra visualization (zoom etc.) and processing
#' (peak picking, integration, phase correction, assignment) from the web browser.
#' 
#' @import hastaLaVista2
#' @export

nmrium <- function(ppm,yr,yi=NULL,roi=range(ppm),cache=file.path(getwd(),"nmrium"),metadata=NULL
                   ,...){
  #wrappers to ease forwarding
  makeServer <- function(...,observeFrequency,dataType,solvent,nucleus,col,meta){
    new("server",...)
  }
  export <- function(...,baseURL,port,path,protocole,rootDir,init){
    hastaLaVista2::exportReIm(...)
  }
  if(!is.null(yr) & !is.null(yi) & !all(dim(yr) == dim(yi))){
    cat(crayon::red("nmr.spectra.processing::nmrium >> ",
                    "the dimensions of yr and yi must match"))
  }
  
  if(!dir.exists(cache)) dir.create(cache)
  fi <- crop(ppm,roi=roi)
  if (!is.null(yr)) yr <- yr[,fi]
  if (!is.null(yi)) yi <- yi[,fi]
  ppm <- ppm[fi]
  
  usenmrium <- export(x=ppm,yr=yr[1,],yi=yi[1,]
                      ,info=as.list(metadata[1,,drop=FALSE]),...)
  n <- nrow(yr)
  for (i in 2:n) 
    usenmrium <- export(x=ppm,yr=yr[i,],yi=yi[i,]
                        ,info=as.list(metadata[i,,drop=FALSE])
                        ,output=usenmrium,...)
  
  fileServer <- makeServer(rootDir=cache,...)
  httpuvServer <- initServer(fileServer, force=TRUE)
  nmriumView(usenmrium, fileServer, openBrowser=TRUE)
  invisible(httpuvServer)
}