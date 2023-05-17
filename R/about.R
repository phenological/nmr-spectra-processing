#' A cheat sheet for ml in js
about.signal <- function() {
  cat("In ml-nmrProcessing, a signal object is the specification of a (potentially) assigned region of a NMR spectrum.
Signals may contain a list of fitte peaks (see about.NMRPeak1D).
Here is the full list of attributes of a signal (all attributes other than delta are optional):
  intensity: number
  shape: Shape1D
  parameters: ParametersFromSignal
  delta: number
  id: string
  js: Jcoupling[]
  atoms: number[]
  assignment: string
  kind: string
  multiplicity: string
  diaIDs: string[]
  nbAtoms: number
  integration: number
  peaks: NMRPeak1D[]
  statistic: {mean: numberm sd: numberm min: number, max: number, nb: number}
On the R side, the object's attributes are specified as named list elements.")
}

about.NMRPeak1D <- function(){
  cat("In ml-nmrProcessing, NMRPeak1D is an object representing a function fit to a peak (singlet) of the spectrum. It has attributes:
  x: number
  y: number
  width: number
  id: string
  kind: string
  shape: Shape1D
On the R side, the object's attributes are specified as named list elements.")
}

about.Shape1D <- function(){
  cat("In ml-nmrProcessing, Shape1D is an alias for any of three objects representing functions that may be used for signal fitting. These are:
      GaussianShape1D: {kind: 'gaussian', fwhm: number}
      LorentizanShape1D: {kind: 'lorentzian', fwhm: number}
      PseudoVoigthShape1D: {kind: 'pseudoVoigt', fwhm: number, mu: number}
     On the R side, the object's attributes are specified as named list elements.")
}

about.JCoupling <- function(){
  cat("In ml-nmrProcessing, JCoupling is an object representing spin coupling. Of all its attributes, only 'coupling' is mandatory:
       coupling: number
       atoms: number[]
       assignment: string | string[]
       diaIDs: string[]
       multiplicity: string
       pathLength: number
      On the R side, the object's attributes are specified as named list elements.")
}