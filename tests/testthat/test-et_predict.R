context("test-et_predict")

test_that("et_predict works", {
  
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600)
  
  freq_range = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'hw95s',
                      freq_range = freq_range)
  et$predict(method = 'gravity')
  expect_equal(et$output$gravity, 
               c(-448.580, -564.521), 
               tolerance = .0001) 
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'ksm03',
                      freq_range = freq_range)
  et$predict(method = 'gravity')
  expect_equal(et$output$gravity, 
               c(-448.648, -564.549), 
               tolerance = .0001)
  
  
})
