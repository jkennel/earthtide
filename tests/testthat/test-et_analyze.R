context("test-et_analyze")

test_that("dimensions correct", {
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600)
  
  wave_groups = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'hw95s',
                      wave_groups = wave_groups)
  
  out <- et$analyze(method = 'gravity')$tide()
  expect_equal(NROW(out), length(tms))
  expect_equal(NCOL(out), 3)
  
  
  tms <- as.POSIXct('1990-01-01', tz = 'UTC')
  wave_groups = data.frame(start = c(0, 1.5), end = c(1.4, 8))
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      wave_groups = wave_groups)
  
  out <- et$analyze(method = 'gravity')$tide()
  expect_equal(nrow(out), length(tms))
  expect_equal(ncol(out), 5)
  
  
  
  
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + 0:(24*31) * 3600
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      wave_groups = wave_groups)
  
  out <- et$analyze(method = 'gravity', astro_update = 24)$tide()
  expect_equal(nrow(out), length(tms))
  expect_equal(ncol(out), 5)
  
  
})
