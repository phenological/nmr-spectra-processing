test_that("shiftSeries works", {
  v <- which.min(abs(ppm-9.8))
  n <- length(ppm)
  savo <- shiftSeries(avo3[1,],30,from=v:n)
  expect_equal(avo3[1,1:(n-30)], savo[31:n])
  expect_lt(max(savo[1:30]),799)
  
  savo <- shiftSeries(avo3[1,],-20,padding="zeroes")
  expect_equal(savo[1:(n-20)], avo3[1,21:n])
  expect_equal(savo[(n-19):n],rep(0,20))
})