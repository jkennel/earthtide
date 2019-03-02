context("delta_t tests")

test_that("delta_t works", {
  expect_equal(delta_t(as.POSIXct('1994-06-01', tz= 'UTC')), 60.3530, 
               tolerance = 0.001)
  
  expect_equal(delta_t(as.POSIXct('1994-07-01', tz= 'UTC')), 60.4012,
               tolerance = 0.001)
})




