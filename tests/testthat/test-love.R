context("love")

test_that("love works", {
  tms <- as.POSIXct('1995-01-01', tz = 'UTC')
  
  wave_groups = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 49.00937,
                      longitude = 8.40444,
                      elevation = 120,
                      gravity = 9.8127, 
                      cutoff = 1e-10,
                      wave_groups = wave_groups)
  
  love_dat <- et$love_params

  # comparison to ETERNA
  expect_equal(13.943036,  love_dat$dom0)
  expect_equal(15.073729,  love_dat$domr)
  expect_equal(-0.000625,  love_dat$dgr)
  expect_equal(-0.002505,  love_dat$dhr)
  expect_equal(-0.001261,  love_dat$dkr)
  expect_equal( 0.0000781, love_dat$dlr)
  expect_equal( 0.001244,  love_dat$dtr)
  
  
})
