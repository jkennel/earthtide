context("test-get_iers")

test_that("get iers works", {
  testthat::skip_on_cran()
  testthat::skip_on_ci()
  testthat::skip_if_offline()

  eop <- get_iers()
  expect_s3_class(eop, "data.frame")
  expect_equal(ncol(eop), 8L)
  expect_gt(nrow(eop), 20990)
})
