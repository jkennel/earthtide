context("test-earthtide")

test_that("earthtide works", {

  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600)
  
  wave_groups = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  et$predict(method = 'tidal_potential', astro_update = 1L)
  et$predict(method = 'tidal_tilt', astro_update = 1L)
  et$predict(method = 'vertical_displacement', astro_update = 1L)
  et$predict(method = 'horizontal_displacement', astro_update = 1L)
  et$predict(method = 'n_s_displacement', astro_update = 1L)
  et$predict(method = 'e_w_displacement', astro_update = 1L)
  et$predict(method = 'vertical_strain', astro_update = 1L)
  et$predict(method = 'areal_strain', astro_update = 1L)
  et$predict(method = 'horizontal_strain', astro_update = 1L)
  et$predict(method = 'volume_strain', astro_update = 1L)
  et$predict(method = 'ocean_tides', astro_update = 1L)
  
  tide <- et$tide()
  expect_equal(tide$tidal_potential, c(1.422, 1.890), tolerance = 0.001)
  expect_equal(tide$tidal_tilt, c(-17.845, -20.578), tolerance = 0.001)
  expect_equal(tide$vertical_displacement, c(87.822, 110.715), tolerance = 0.001)
  expect_equal(tide$horizontal_displacement, c(-47.34171, -55.21554), tolerance = 2) # from solidearthtide
  expect_equal(tide$n_s_displacement, c(-47.34171, -55.21554), tolerance = 2) # from solidearthtide
  expect_equal(tide$e_w_displacement, c(31.11908, 13.03227), tolerance = 2) # from solidearthtide
  expect_equal(tide$vertical_strain, c(-5.533, -6.964), tolerance = 0.001)
  expect_equal(tide$areal_strain,  c(16.600, 20.893), tolerance = 0.001)
  expect_equal(tide$volume_strain, c(11.067, 13.929), tolerance = 0.001)
  expect_equal(tide$horizontal_strain, c(7.958, 11.823), tolerance = 0.001)
  expect_equal(tide$ocean_tides, c(146.327, 183.187), tolerance = 0.01)
  
  # expect_equal(et$predict(method = 'horizontal_displacement', astro_update = 1L)[1:2],
  #              c(87.822, 110.715),
  #              tolerance = .001)
  
  # non-equal spacing
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600, 7000)
  wave_groups = data.frame(start = 0, end = 8)
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)

  expect_warning(et$predict(method = 'tidal_potential', astro_update = 2L))
  expect_silent(et$predict(method = 'tidal_potential', astro_update = 1L))
  

  expect_silent(et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-4,
                      catalog = 'ksm04',
                      wave_groups = wave_groups, 
                      eop = earthtide:::dut1))
  
  expect_error(et <- Earthtide$new(utc = tms, 
                              latitude = 52.3868,
                              longitude = 9.7144,
                              elevation = 110,
                              gravity = 9.8127, 
                              cutoff = 1.0e-4,
                              catalog = 'ksm04',
                              wave_groups = wave_groups, 
                              eop = data.frame(datetime = 0)))
  
  
  mult <- 1.5
  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + c(0, 3600)
  wave_groups = data.frame(start = 0, end = 8, multiplier = 1.5)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  et$predict(method = 'tidal_potential', astro_update = 1L)
  
  tide <- et$tide()
  expect_equal(tide$tidal_potential, 1.5 * c(1.422, 1.890), tolerance = 0.001)
  
  
  
  
  wave_groups <- na.omit(eterna_wavegroups[eterna_wavegroups$time == 'all', c('start', 'end')])
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-5,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  
  wave_groups <- na.omit(eterna_wavegroups[eterna_wavegroups$time == 'all', c('start', 'end')])
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-5,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  et$predict(method = 'tidal_potential', astro_update = 1L)
  et$lod_tide()
  et$pole_tide()
  et_r6 <- et$tide()
  
  et_fun <- calc_earthtide(utc = tms,
                      do_predict = TRUE,
                      method = c('tidal_potential', 'lod_tide', 'pole_tide'),
                      astro_update = 1,
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-5,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  expect_equal(et_r6, et_fun)

  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-5,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  et_mat <- et$predict(method = 'tidal_potential', astro_update = 1L,
                       return_matrix = TRUE)
  expect_equivalent(as.matrix(et_r6[, c('tidal_potential')]), et_mat)
  
  et_mat2 <- calc_earthtide(utc = tms,
                           do_predict = TRUE,
                           method = c('tidal_potential'),
                           astro_update = 1,
                           latitude = 52.3868,
                           longitude = 9.7144,
                           elevation = 110,
                           gravity = 9.8127, 
                           cutoff = 1.0e-5,
                           catalog = 'ksm04',
                           wave_groups = wave_groups,
                           return_matrix = TRUE)
  
  expect_equivalent(et_mat2, et_mat)
  
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-5,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  et$analyze(method = 'tidal_potential', astro_update = 1L)
  et_r6 <- et$tide()
  
  et_fun <- calc_earthtide(utc = tms,
                      do_predict = FALSE,
                      method = c('tidal_potential'),
                      astro_update = 1,
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-5,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  expect_equal(et_r6, et_fun)
  
  
  
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 52.3868,
                      longitude = 9.7144,
                      elevation = 110,
                      gravity = 9.8127, 
                      cutoff = 1.0e-5,
                      catalog = 'ksm04',
                      wave_groups = wave_groups)
  
  et_mat <- et$analyze(method = 'tidal_potential', astro_update = 1L,
                       return_matrix = TRUE)
  expect_equivalent(as.matrix(et_r6[, -c(1)]), et_mat)
  
  et_mat2 <- calc_earthtide(utc = tms,
                 do_predict = FALSE,
                 method = c('tidal_potential'),
                 astro_update = 1,
                 latitude = 52.3868,
                 longitude = 9.7144,
                 elevation = 110,
                 gravity = 9.8127, 
                 cutoff = 1.0e-5,
                 catalog = 'ksm04',
                 wave_groups = wave_groups,
                 return_matrix = TRUE)
  expect_equivalent(et_mat2, et_mat)
  
  

  
  
  expect_error(calc_earthtide(utc = tms,
                           do_predict = FALSE,
                           method = c('tidal_potential', 'gravity'),
                           astro_update = 1,
                           latitude = 52.3868,
                           longitude = 9.7144,
                           elevation = 110,
                           gravity = 9.8127, 
                           cutoff = 1.0e-5,
                           catalog = 'ksm04',
                           wave_groups = wave_groups))
  
  
  
  

  tms <- as.POSIXct('1990-01-01', tz = 'UTC') + seq(0, 86400*5, 3600)
  
  
  et_mat2 <- calc_earthtide(utc = tms,
                            do_predict = FALSE,
                            method = c('tidal_potential'),
                            astro_update = 1,
                            latitude = 52.3868,
                            longitude = 9.7144,
                            elevation = 110,
                            gravity = 9.8127, 
                            cutoff = 1.0e-5,
                            catalog = 'ksm04',
                            wave_groups = wave_groups,
                            return_matrix = TRUE)

  et_mat3 <- calc_earthtide(utc = tms,
                            do_predict = FALSE,
                            method = c('tidal_potential'),
                            astro_update = 10,
                            latitude = 52.3868,
                            longitude = 9.7144,
                            elevation = 110,
                            gravity = 9.8127, 
                            cutoff = 1.0e-5,
                            catalog = 'ksm04',
                            wave_groups = wave_groups,
                            return_matrix = TRUE)
  expect_equivalent(et_mat2, et_mat3)
  
  
  
})

