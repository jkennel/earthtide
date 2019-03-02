context("test-print")

test_that("print works", {

  wave_groups = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = as.POSIXct('2000-01-01', tz = 'UTC'), 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'hw95s',
                      wave_groups = wave_groups)
  
  expect_output(et$print(), regexp = NULL)

  
})
