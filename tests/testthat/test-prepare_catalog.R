context("test-prepare_catalog")

test_that("prepare_catalog works", {
  
  expect_silent(.prepare_catalog(1e-10, wave_groups = NULL, catalog = 'ksm04'))
  
  wave_groups <- data.frame(start = 1, end = 2, multiplier = 1.0)
  expect_silent(.prepare_catalog(1e-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- data.frame(start = 1, end = 2)
  expect_silent(.prepare_catalog(1e-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- data.frame(a = 1, b = 2)
  expect_silent(.prepare_catalog(1e-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- data.frame(a = 1, b = 2, c = 1.0)
  expect_silent(.prepare_catalog(1e-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- data.frame(a = c(1,3), b = c(2,4), c = c(1.0, 1.0))
  expect_silent(.prepare_catalog(1e-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- data.frame(a = c(1,3), b = c(2,4), c = c(1.0, 1.0), d = c(2, 5))
  expect_warning(.prepare_catalog(1-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- data.frame(start = 1)
  expect_error(.prepare_catalog(1e-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- na.omit(eterna_wavegroups[eterna_wavegroups$time == 'all', c('start', 'end')])
  expect_silent(.prepare_catalog(1e-10, wave_groups = wave_groups, catalog = 'ksm04'))
  
  wave_groups <- na.omit(eterna_wavegroups[eterna_wavegroups$time == 'all', c('start', 'end')])
  expect_silent(.prepare_catalog(1e-4, wave_groups = wave_groups, catalog = 'ksm04'))
  
  
})
