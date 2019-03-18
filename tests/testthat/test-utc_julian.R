context("test-utc_julian")

test_that("multiplication works", {
  expect_equal(utc_julian(as.POSIXct('2019-01-01 00:00:00', tz = 'UTC')), 2458484.50000)
  expect_equal(utc_julian(as.POSIXct('2000-01-01 12:00:00', tz = 'UTC')), 2451545.000000)
  
  expect_equal(mod_julian_utc(utc_mod_julian(as.POSIXct('2000-01-01 12:00:00', tz = 'UTC'))), 
               as.POSIXct('2000-01-01 12:00:00', tz = 'UTC'))
  
})
