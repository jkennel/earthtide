context("test-et_analyze")

test_that("dimensions correct", {
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600)
  
  freq_range1 = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'hw95s',
                      freq_range = freq_range1)
  
  out1 <- et$analyze(method = 'gravity')
  expect_equal(NROW(out1), length(tms))
  expect_equal(NCOL(out1), 3)
  
  
  tms <- as.POSIXct('1990-01-01', tz = 'UTC')
  freq_range2 = data.frame(start = c(0, 1.5), end = c(1.4, 8))
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'ksm03',
                      freq_range = freq_range2)
  
  out2 <- et$analyze(method = 'gravity')
  expect_equal(nrow(out2), length(tms))
  expect_equal(ncol(out2), 5)
  
})
