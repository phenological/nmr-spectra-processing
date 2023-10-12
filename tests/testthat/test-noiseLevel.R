test_that("noiseLevel works", {
  expect_equal(round(noiseLevel(ppm, avo3[1,]),4), 540.5311)
  expect_equal(round(noiseLevel(ppm, avo3, level=0.95, rOref=c(-0.1,-0.03)),4)
               ,c(488.2212, 325.0211, 425.8497))
})