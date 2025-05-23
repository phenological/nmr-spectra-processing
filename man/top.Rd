% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{top}
\alias{top}
\title{Get the \emph{n} spectra with the highest intensity on the given chemical shift
range or value}
\usage{
top(ppm, Y, cshift, n = 10L, roi = c(-Inf, Inf), bottom = FALSE, index = FALSE)
}
\arguments{
\item{ppm, }{numeric, chemical shift scale.}

\item{Y, }{matrix, numeric, intensities, spectra in rows.}

\item{cshift, }{numeric, optional, query chemical shift value.}

\item{n, }{integer, number of spectra to be returned, 10 by default.}

\item{roi, }{numeric, optional, length 2 vector with the limits of the query
chemical shift range.}

\item{bottom, }{logic, if TRUE returns the \code{n} spectra with the \emph{lowest}
intensity instead.}

\item{index, }{logic, determnies whether to return the row indices of top
spectra (TRUE) or the spectra themselves (FALSE, default).}
}
\value{
If index is TRUE, returns a vector with the row indices of the
\emph{n} spectra with the highest (lowest if \code{bottom}) intensity at
the given \code{cshift} or  within the given \code{roi}. Otherwise (default)
returns a matrix with the corresponding spectra.
}
\description{
Get the \emph{n} spectra with the highest intensity on the given chemical shift
range or value
}
\details{
If both a precise chemical shift and a chemical shift range are
passed, \code{cshift} takes priority
}
