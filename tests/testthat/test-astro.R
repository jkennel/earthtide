context("test-astro")

test_that("astro works", {
  
  tms <- as.POSIXct('1995-01-01', tz = 'UTC')
  
  wave_groups = data.frame(start = 0, end = 8)
  
  et <- Earthtide$new(utc = tms, 
                      latitude = 49.00937,
                      longitude = 8.40444,
                      elevation = 120,
                      gravity = 9.8127, 
                      cutoff = 1.0e-10,
                      wave_groups = wave_groups)
  
  ast <- 
    c(16.94984832601,
      271.63781389767, 
      280.18224962285,
      239.87476355745,
      138.23542575910,
      282.85135750349,
      337.52710141118,
      135.61974634875,
      118.21753519044,
      242.51565706409,
      348.89354141160)
  
  ast_der <- 
    c(14.49205211998,
       0.54901651990,
       0.04106863988,
       0.00464181459,
       0.00220640711,
       0.00000196146,
       0.17051571086,
       0.06675703048,
       0.02183629516,
       0.00346372662,
       0.00139574608)
  
  expect_equal(as.numeric(et$astro$astro), ast, 
               tolerance = 1e-4)
  expect_equal(as.numeric(et$astro$astro_der), ast_der,
               tolerance = 1e-4)
  
})
