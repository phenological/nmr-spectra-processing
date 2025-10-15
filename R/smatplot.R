#' Fast multiple spectra overlay plot
#' 
#' Wrapper to \code{\link[graphics]{matplot}} for plotting spectra. Crude but 
#' fast compared to fancier alternatives. Allows to split spectra in multiple
#' plots, adjust plot resolution to screen resolution, choose between two layouts,
#' add admittedly crude legends, and (nearly) full \code{\link[graphics]{graphics}}
#' customization
#' 
#' @param ppm numeric, spectra chemical shift scale a.k.a ppm
#' @param y numeric vector or matrix, nmr intensities, spectra in rows. If neither
#' vector nor matrix, argument is cast as.matrix. Type compliance is recommended.
#' @param roi numeric, optional, limits of the Region of Interest to be plotted.
#'  Defaults to the range of \code{ppm}.
#' @param by numeric, optional, number of spectra to be overlayed on each plot
#' @param lty optional, line type, defaults to 1 (continuous line). Note that this
#' parameter is passed both to \code{\link[graphics]{matplot}} and \code{\link[graphics]{legend}},
#' so you always get lines in your legend (see ... below)
#' @param legend, optional, position of the legend as specified in
#' \code{\link[graphics]{legend}}. No legend is printed by default.
#' @param reverse logic, optional, if TRUE (default) the x scale of the plot 
#' increases from right to left as it's customary in NMR spectroscopy
#' @param resolution, character. If "full" all data are passed to 
#' \code{\link[graphics]{matplot}}. If "dev", data are binned to fit the graphic
#'  device's resolution before plotting. By default \code{smatplot} chooses
#'  according to data size. See Details.
#' @param reduce, function used to compute bin values. See Details
#' @param col vector to color by. Works just as in R base plotting. The default is \emph{Set1} copied from 
#' \code{\link[RColorBrewer]{RColorBrewer}}.
#' @param stacked logical, if TRUE plots spectra stacked along the y axis.
#' If FALSE (default) overlays spectra on top of each other.
#' @param ..., additional arguments for customization. These arguments are
#' passed either to \code{\link[graphics]{matplot}} for plot customization, or
#' to  \code{\link[graphics]{legend}} for legend customization. Beware that some
#' of these parameters are shared, such as \code{ltw} and the non-optional
#' parameters \code{lty} and \code{col}. This means, for instance, that it is 
#' not possible to get a legend without lines or to unmatch legend linewidths 
#' from the plot's.
#' @returns NULL
#' @details
#' tl;dr: For spectra with less than ~1.2 million points the default method
#' renders fast and exact plots. Otherwise:
#' * the default method adds a quick pre-processing step to make the plot fast 
#' to render, at the cost of inaccurate peak heights
#' * \code{resolution="dev", reduce=max} pre-processes slower, renders as fast 
#' and gives accurate peak heights, but introduces ppm inaccuracies visible when 
#' zooming in
#' * \code{resolution="full"} is exact but very slow to render.
#' 
#' Further details follow.
#' 
#' When working with high resolution spectra it is likely that the resolution
#' of the spectrum is higher than the pixel resolution of the graphic device,
#' which places an unnecessary burden on the renderer. By default \code{smatplot}
#' bins the spectra matrix to match the pixel resolution of the active graphic 
#' device (as reported by \code{\link[grDevices]{dev.size}}) if the input spectra
#' contain more than 1.2 million points total. The user can force full resolution
#' or binned resolution by setting \code{resolution} to "full" or "dev" 
#' respectively.
#' 
#' By default, bin values are computed by sampling each spectrum at regular 
#' intervals. This method is fast, which is good for the visualization of large 
#' datasets, but it produces inaccurate estimates of peak height. Alternatively, 
#' the user may pass a function to the \code{reduce} parameter that will be used 
#' to compute the bins. In this case, each spectrum is partitioned with 
#' \code{\link[base]{split}} and the \code{reduce} function is applied to compute 
#' each bin's intensity value and ppm. This adds a expensive pre-processing
#' step that makes the plot slower to compute but the result renders much 
#' faster than \code{resolution="full"}. Beware that spectra rendered in this manner
#' are chemical shift-accurate only up to the current graphic resolution.
#' This means that inaccuracies in the frequency coordinate may appear
#' if you enhance the resolution of the graphic device after plotting e.g. by
#' enlarging the plot window in RStudio.
#' From experience,\code{reduce=max} produces the most accurate picture of peak 
#' intensities, while \code{reduce=mean} or \code{median} provide a compromise 
#' between intensity and chemical shift accuracy.
#' 
#' @importFrom graphics matplot
#' @importFrom grDevices dev.size
#' @export
smatplot <- function (ppm, y, roi, by, lty = 1, legend, reverse = TRUE, 
  resolution = c("full", "dev", ifelse(length(as.matrix(y)) > 
    1200000, "dev", "full"))[3], reduce, col = c("#E41A1C", 
    "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", 
    "#A65628", "#F781BF", "#999999"), stacked = FALSE, ...) 
{
  lmatplot <- function(..., fill, border, angle, density, 
    bty, bg, box.lwd, box.lty, box.col, pt.bg, pt.cex, pt.lwd, 
    xjust, yjust, x.intersp, y.intersp, adj, text.width, 
    text.col, text.font, merge, trace, plot, ncol, horiz, 
    title, insert, xpd, title.col, title.adj, title.cex, 
    title.font, seg.len) {
    matplot(...)
  }
  llegend <- function(..., adj, ann, cex.axis, cex.lab, cex.main, 
    cex.sub, col.axis, col.lab, col.main, col.sub, crt, 
    err, family, fg, font, font.axis, font.lab, front.main, 
    font.sub, lab, las, lend, ljoin, lmitre, mai, mar, mex, 
    mgp, mkh, page, smo, srt, tck, xaxp, xaxs, xaxt, xpd, 
    yaxp, yaxs, yaxt, type, main, sub, xlab, ylab, asp, 
    xlim, ylim, log, add, verbose, frame) {
    legend(...)
  }
  
  
  if (!is.matrix(y) & !is.vector(y)) {
    cat(crayon::yellow("nmr-spectra-processing::smatplot >>", 
      "Argument Y being cast as.matrix.\n", "Unpredictable results may follow if casting to", 
      "numeric matrix fails\n"))
    y <- as.matrix(y)
  }
  if (is.matrix(y)) 
    y <- t(y)
  else y <- as.matrix(y)
  
  #Need to cycle colors myself to achieve consistency when using the "by" parameter
  n <- ncol(y) %/% length(col)
  k <- ncol(y) %% length(col)
  col <- c(rep(col,n),col[1:k])
  
  if (!missing(roi)) {
    roi <- sort(roi)
    fi <- ppm >= roi[1] & ppm <= roi[2]
    ppm <- ppm[fi]
    y <- y[fi, , drop = FALSE]
  }
  else {
    roi <- range(ppm)
  }
  pointz <- dim(y)[1L]
  if (resolution == "dev") {
    pixels <- grDevices::dev.size(units = "px")[1]
    if (pointz > pixels) {
      cat(crayon::yellow("nmr.spectra.processing::smatplot >>", 
        "Binning spectra to fit the graphic device's resolution\n", 
        "If you want to keep full resolution set argument resolution='full'\n"))
      pointsPerPixel <- pointz%/%pixels
      if (missing(reduce)) {
        fi <- 1:pointz%%pointsPerPixel == 0
        ppm <- ppm[fi]
        y <- y[fi, , drop = FALSE]
      }
      else {
        cat(crayon::yellow("nmr.spectra.processing::smatplot >>", 
          "Applying your custom reduce function\n"))
        bins <- as.factor(1:pointz%/%pointsPerPixel)
        ppm <- unname(sapply(split(ppm, bins), reduce))
        y <- apply(y, 2, function(v) {
          unname(sapply(split(v, 1:pointz%/%pointsPerPixel), 
            reduce))
        })
      }
    }
  }
  if (reverse) {
    roi <- rev(roi)
  }
  extrargs <- list(...)
  if (!"xlab" %in% names(extrargs)) 
    extrargs$xlab <- "ppm"
  if (!"ylab" %in% names(extrargs)) 
    extrargs$ylab <- NA
  n <- ncol(y)
  if (missing(by)) {
    groups <- list(1:n)
  }
  else {
    groups <- split(1:n, ceiling(1:n/by))
  }
  if ("ylim" %in% names(extrargs)) 
    vwindow <- diff(extrargs$ylim)
  for (g in groups) {
    yg <- y[, g, drop = FALSE]
    colg <- col[g]
    if (stacked) {
      if ("ylim" %in% names(extrargs)) {
        extrargs$ylim <- c(extrargs$ylim[1], extrargs$ylim[1] + 
          vwindow * length(g))
      }
      else {
        vwindow <- max(yg)
      }
      for (i in 1:length(g)) yg[, i] <- yg[, i] + vwindow * 
        (i - 1)
    }
    do.call(lmatplot, c(list(x = ppm, y = yg, type = "l", 
      xlim = roi, lty = lty, col = colg), extrargs))
    if (!missing(legend)) {
      if(is.factor(col)){
        do.call(llegend, c(list(x = legend, legend = levels(col), 
          col = 1:length(levels(col)), lty = lty), extrargs))
      }
      else {
        do.call(llegend, c(list(x = legend, legend = unique(col), 
          col = unique(col), lty = lty), extrargs))
      }
    }
  }
}