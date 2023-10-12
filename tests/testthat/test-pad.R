test_that("pad works", {
  av <- avo3[1,]
  pav <- pad(av,7,side=-1)
  expect_equal(length(pav),40813L)
  expect_equal(pav[1:40806],avo3RandShift[1,])
  
  pav <- pad(av,8,side=0,method="sampling",from=39933:40806)
  expect_equal(length(pav),40806+16)
  expect_equal(pav[9:40814],avo3[1,])
  expect_lt(max(pav[c(1:8,40815:40822)]),799)
  
  pav <- pad(av,10,side=1,method="circular")
  expect_equal(pav[40807:40816],av[1:10])
  pav <- pad(av,10,side=-1,method="circular")
  expect_equal(pav[1:10],av[40797:40806])
})