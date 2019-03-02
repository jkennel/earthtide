context("test-utc_julian_2000")

test_that("utc_julian_2000 works", {
  expect_equal(utc_julian_2000(as.POSIXct('2000-01-01 12:00:00', tz = 'UTC')), 0)
  expect_equal(utc_julian_2000(as.POSIXct('2019-01-01 12:00:00', tz = 'UTC')), 599616000/86400/36525 )
})
