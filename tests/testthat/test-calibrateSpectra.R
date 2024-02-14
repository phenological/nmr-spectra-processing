test_that("calibration works", {
  calib <- calibrateSpectra(ppm,avo3RandShift,"tsp")
  fi <- crop(ppm,-0.02,0.02)
  expect_equal(calib[,fi],avo3[,fi])

    expect_equal(calibrateSpectra(ppm,avo3,"alanine",frequency=400), avo3CalAla)
    
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
    calib <- calibrateSpectra(ppm, avo3, sucrose
                              ,frequency=400, rOref=c(5.4,5.435))
    expect_equal(calib,avo3CalSucro)
    
    sucrosecal <- calibrateSignal(ppm,avo3[1,],sucrose)
    expect_equal(round(sucrosecal@chemicalShift,5),5.41577)
    expect_equal(round(sapply(sucrosecal@peaks,function(p) p@x),5)
                 ,c(5.41107,5.42047))
    expect_equal(round(sapply(sucrosecal@peaks,function(p) p@y),5)
                 ,c(1L,1L))
    
    expect_no_condition(calibrateSpectra(ppm,avo3[1,]))
    expect_equal(calibrateSpectra(as.character(ppm),avo3,"alanine",frequency=400)
                 ,avo3CalAla)
    expect_equal(calibrateSpectra(ppm,as.character(avo3[1,]),"alanine",frequency=400)
    ,avo3CalAla[1,])
})

