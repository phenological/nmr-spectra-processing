test_that("signalToY works", {
  sucrose <- new("NMRSignal1D"
                 ,id="Sucrobation"
                 ,chemicalShift=5.416
                 ,shape=list(name="pseudoVoigt"
                             ,params = list(mu=0.85,fwhm=0.0015,base=0)
                 )
                 ,peaks = list(new("NMRPeak1D",x=5.4113, y=1)
                               ,new("NMRPeak1D",x=5.4207, y=1)
                 )
  )
  handmade <- pseudoVoigt(ppm,5.4113,1
                          ,0.0015,0.85) + pseudoVoigt(ppm,5.4207,1,0.0015,0.85)
  expect_equal(signalToY(sucrose,ppm),handmade)
  
})