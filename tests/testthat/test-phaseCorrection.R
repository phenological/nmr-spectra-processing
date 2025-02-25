test_that("phaseCorrection works", {
  
  x <- (0:128) * 2 * pi / 128
  yr <- cos(x)
  yi <- sin(x)
  
  res <- phaseCorrection(yr, yi, 0, 0)
  
  expect_equal(yr, res[[1]])
  expect_equal(yi, res[[2]])
  
  res <- phaseCorrection(yr, yi, 180, 0)
  
  expect_equal(yr, -res[[1]])
  expect_equal(yi, -res[[2]])
  
  res <- phaseCorrection(yr, yi, 360, 0)
  
  expect_equal(yr, res[[1]])
  expect_equal(yi, res[[2]])
  
  res <- phaseCorrection(yr, yi, -90, 90)
  
  x <- (0:128) * (6.0/4.0) * pi / 128.0 
  yr_ex <- cos(x)
  yi_ex <- sin(x)

  expect_equal(yr_ex, res[[1]], tolerance=1e-2)
  expect_equal(yi_ex, res[[2]], tolerance=1e-2)

})