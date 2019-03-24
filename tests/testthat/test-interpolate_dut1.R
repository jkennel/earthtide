context("interpolate_dut1 tests")

test_that("interpolate_dut1 works", {
  
  expect_equal(interpolate_dut1(as.POSIXct('1994-06-01', tz= 'UTC'), 'ddt', dut1),
               60.3530, tolerance = 0.001)
  expect_equal(interpolate_dut1(as.POSIXct('1994-07-01', tz= 'UTC'), 'ddt', dut1), 
               60.4012, tolerance = 0.001)
  
  
  expect_equal(interpolate_dut1(as.POSIXct('2018-12-02', tz = 'UTC'), 'ut1_utc', dut1), 
               -0.0137963, tolerance = 1e-3)
  expect_equal(interpolate_dut1(as.POSIXct('2019-01-01', tz = 'UTC'), 'ut1_utc', dut1),
               -0.0361691, tolerance = 1e-3)
  
  
  
})




