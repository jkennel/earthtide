context("test-prepare_station")

test_that("prepare_station works", {
  tms <- as.POSIXct('1995-01-01', tz = 'UTC')
  
  wave_groups = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 49.00937,
                      longitude = 8.40444,
                      elevation = 120,
                      gravity = 9.8127, 
                      cutoff = 1e-10,
                      wave_groups = wave_groups)
  
  geodetic <- et$station$dgz
  eterna_result <- c(
    246.45038826,
    600.43570011,
    261.76537479,
    -75.76801831,
    917.50564499,
    782.14902410,
    278.72506065,
    -644.30912324,
    711.47974781,
    1344.49043258,
    838.93750710,
    259.05392802,
    -1077.54405654, 
    -152.58086909,
    1511.35174586,
    1577.22960769,
    807.87471610,
    223.20077736,
    -971.53577214,
    -1335.91779371,
    924.46283635,
    2113.43269097,
    1636.58158318,
    726.27093065,
    183.21337179
  )
  expect_equal(geodetic,  eterna_result)
  
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 49.00937,
                      longitude = 8.40444,
                      elevation = 120,
                      gravity = 0, 
                      cutoff = 1e-10,
                      wave_groups = wave_groups)
  grav <- et$station$gravity
  expect_equal(grav,  gravity_station(49.00937, 120))
})
