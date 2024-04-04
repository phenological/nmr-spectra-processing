test_that("normalize works", {
  expect_equal(unique(rowSums(normalize(avo3,sum))),1)
  expect_equal(unique(max(normalize(avo3,max))),1)
  expect_equal(round(pqn(avo3),3),c(1.000,1.017,0.896))
})