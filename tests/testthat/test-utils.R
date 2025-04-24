test_that("crop works", {
  expect_equal(avo3[,crop(ppm,-0.02,0.02)], avo3tsp)
  expect_equal(avo3[,crop(ppm,roi=c(1.47,1.52))], avo3Ala)
  expect_equal(avo3[,crop(ppm,end=1.52)], avo3StartToAla)
  expect_equal(avo3[,crop(ppm,start=5.4)], avo3SucrToEnd)
})
test_that("getI works", {
  expect_equal(getI(ppm,5),19422L)
  expect_equal(getI(ppm,c(5,5.2))
               ,c(19422L,19858L))
})

test_that("top works", {
  expect_equal(top(ppm, avo3, 1.5, n=1, index=TRUE), 3L)
  expect_equal(top(ppm, avo3, 1.5), avo3[3:1,])
  expect_equal(top(ppm,avo3,roi=c(1.503,1.51),n=1, index=TRUE), 2L)
})

# test_that("normalizeSignal works", {
#   sucrose <- new("NMRSignal1D"
#                  ,id="Sucro"
#                  ,chemicalShift=5.416
#                  ,shape=list(name="pseudoVoigt"
#                              ,params = list(mu=0.85,fwhm=0.0015,base=0)
#                  )
#                  ,peaks = list(new("NMRPeak1D",x=5.4113, y=123)
#                                ,new("NMRPeak1D",x=5.4207, y=125)
#                  )
#   )
#   expect_s4_class(normalizeSignal(sucrose),"NMRSignal1D")
#   summit <- max(sapply(normalizeSignal(sucrose)@peaks,function(p) p@y))
#   summitI <- which.max(sapply(normalizeSignal(sucrose)@peaks,function(p) p@y))
#   expect_equal(summit,1)
#   expect_equal(summitI,2)
# })