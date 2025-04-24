#naive unidecnmr test, using pseudoVoigts
unidec <- function(ppm,I,fwhm=1/600,roi
                   ,noise = noiseLevel(ppm,I,rOref = c(9.9,10))
                   ,maxIt = 10){
  # noise <- noise * 
  fi <- crop(ppm,roi=roi)
  ppm <- ppm[fi]
  I <- I[fi]
  n <- length(I)
  f <- I
  step <- 0
  while (step < maxIt){
    step <- step + 1
    peaks <- sapply(1:n, function(i) fpseudoVoigt(ppm,mu=0.85,mean=ppm[i],max=f[i],fwhm=fwhm))
    reconstruction <- rowSums(peaks)
    f <- f * I / reconstruction
  }
  smatplot(ppm,reconstruction,roi=roi)
  smatplot(ppm,I,roi=roi,add=T,col="black")
  apply(peaks,2,function(p) lines(ppm,p,type="h",col="blue"))
}
