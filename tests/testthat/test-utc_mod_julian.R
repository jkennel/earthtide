context("test-utc_mod_julian")

test_that("utc_mod_julian works", {
  expect_equal(utc_mod_julian(as.POSIXct('1995-10-10 00:00:00', tz = 'UTC')), 50000)
  expect_equal(utc_mod_julian(as.POSIXct('2000-01-01 12:00:00', tz = 'UTC')), 51544.5)
  
  expect_equal(
    julian_mod_julian(
      mod_julian_julian(
        utc_mod_julian(as.POSIXct('2000-01-01 12:00:00', tz = 'UTC'))
      )
    ),
    51544.5)
  
  
})
