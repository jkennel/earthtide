context("test-ut1_utc")

test_that("ut1_utc works", {
  
  expect_equal(ut1_utc(as.POSIXct('2018-12-02', tz = 'UTC')), -0.0137963, tolerance = 1e-3)
  expect_equal(ut1_utc(as.POSIXct('2019-01-01', tz = 'UTC')), -0.0361691, tolerance = 1e-3)

})
