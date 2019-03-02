context("test-earthtide")

test_that("earthtide works", {
  
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600)
  
  freq_range = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'ksm03',
                      freq_range = freq_range)
  
  et$predict(method = 'tidal_potential', astro_update = 1L)
  et$predict(method = 'tidal_tilt', astro_update = 1L)
  et$predict(method = 'vertical_displacement', astro_update = 1L)
  et$predict(method = 'vertical_strain', astro_update = 1L)
  et$predict(method = 'areal_strain', astro_update = 1L)
  et$predict(method = 'volume_strain', astro_update = 1L)
  et$predict(method = 'ocean_tides', astro_update = 1L)
  
  expect_equal(et$output$tidal_potential, c(1.422, 1.890), tolerance = 0.001)
  expect_equal(et$output$tidal_tilt, c(-17.845, -20.578), tolerance = 0.001)
  expect_equal(et$output$vertical_displacement, c(87.822, 110.715), tolerance = 0.001)
  expect_equal(et$output$vertical_strain, c(-5.533, -6.964), tolerance = 0.001)
  expect_equal(et$output$areal_strain,  c(16.600, 20.893), tolerance = 0.001)
  expect_equal(et$output$volume_strain, c(11.067, 13.929), tolerance = 0.001)
  expect_equal(et$output$ocean_tides, c(146.327, 183.187), tolerance = 0.01)
  
  # expect_equal(et$predict(method = 'horizontal_displacement', astro_update = 1L)[1:2],
  #              c(87.822, 110.715),
  #              tolerance = .001)
  
  # non-equal spacing
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600, 7000)
  freq_range = data.frame(start = 0, end = 8)
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'ksm03',
                      freq_range = freq_range)

  expect_warning(et$predict(method = 'tidal_potential', astro_update = 1L))
  
})

