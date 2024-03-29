% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/baselineCorrection.R
\name{baselineCorrection}
\alias{baselineCorrection}
\title{Baseline Correction}
\usage{
baselineCorrection(y, ...)
}
\arguments{
\item{y}{numeric or matrix, single spectrum intensities or intensities matrix
with spectra in rows}

\item{...}{other arguments to \code{\link[ptw]{asysm}}}
}
\value{
a vector or matrix of the same dimensions as the input with
baseline-corrected spectra.
}
\description{
Wrapper to \code{\link[ptw]{baseline.corr}}. Estimates baseline by asymmetric
least squares and subtracts it from the spectra
}
\details{
Presently this is \code{\link[ptw]{baseline.corr}} straight out of
the box, wrapped and re-packaged for convenience.
}
