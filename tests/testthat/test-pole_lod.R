context("test-pole_lod")

test_that("pole and LOD tides work", {
  
  tms <- as.POSIXct('2011-01-01', tz = 'UTC')
  wave_groups = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 14.8861,
                      longitude = 103.5203,
                      elevation = 20,
                      gravity = 0, 
                      cutoff = 1.0e-10,
                      catalog = 'hw95s',
                      wave_groups = wave_groups)
  et$pole_tide()
  et$lod_tide()

  tide <- et$tide()
  
  # This needs to be checked
  # values from /media/kennel/Data/tmp/ET34-X-V71/ET34-ANA-V71/PROJECTS/TEST/TEST.dat
  expect_equal(tide$pole_tide, -21.651, tolerance = 0.05)
  expect_equal(tide$lod_tide, 0.248,  tolerance = 0.05)
  # plot(et$pole_tide(), type='l')
  # plot(et$lod_tide(), type='l')
})
