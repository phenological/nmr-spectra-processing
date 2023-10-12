test_that("integral works", {
  sucrose <- new("NMRSignal1D"
                 ,id="Sucrobation"
                 ,chemicalShift=5.416
                 ,shape=list(name="pseudoVoigt"
                             ,params = list(mu=0.85,fwhm=0.0015,base=0)
                 )
                 ,peaks = list(new("NMRPeak1D",x=5.4113, y=4200)
                               ,new("NMRPeak1D",x=5.4207, y=4200)
                               )
                 )
  mu <- .85
  fwhm <- .0015
  h <- 4200
  expect_equal(integral(sucrose,"fwhm")
               ,2 * mu * fwhm * h * pi/2 + 2 * (1-mu) * fwhm * h * 1.064467
  )
  expect_equal(round(integral(sucrose,"rect",offset=10),5),18.40495)
  expect_equal(integral(list(x=ppm,y=avo3[1,]),"rect")
               , sum(avo3[1,]) * (ppm[2]-ppm[1])
               )
})