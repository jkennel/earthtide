context("test-et_predict")

test_that("et_predict works", {
  tms <- as.POSIXct("1990-01-01", tz = "UTC") + c(0, 3600)

  wave_groups <- data.frame(start = 0, end = 8)

  et <- Earthtide$new(
    utc = tms,
    latitude = 52.3868,
    longitude = 9.7144,
    elevation = 110,
    gravity = 9.8127,
    cutoff = 1.0e-10,
    catalog = "hw95s",
    wave_groups = wave_groups
  )
  tmp <- et$predict(method = "gravity")

  tide <- et$tide()



  expect_equal(tide$gravity,
    c(-448.580, -564.521),
    tolerance = .0001
  )

  et <- Earthtide$new(
    utc = tms,
    latitude = 52.3868,
    longitude = 9.7144,
    elevation = 110,
    gravity = 9.8127,
    cutoff = 1.0e-10,
    wave_groups = wave_groups
  )
  et$predict(method = "gravity")
  tide <- et$tide()

  expect_equal(tide$gravity,
    c(-448.648, -564.549),
    tolerance = .0001
  )

  tms <- as.POSIXct("1990-01-01", tz = "UTC") + 0:(24 * 31) * 3600
  et <- Earthtide$new(
    utc = tms,
    latitude = 52.3868,
    longitude = 9.7144,
    elevation = 110,
    gravity = 9.8127,
    cutoff = 1.0e-10,
    wave_groups = wave_groups
  )

  et$predict(method = "gravity")
  tide <- et$tide()

  expect_equal(tide$gravity[1:2],
    c(-448.648, -564.549),
    tolerance = 0.01
  )


  tms <- as.POSIXct("1990-01-01", tz = "UTC") + 0:1800
  et <- Earthtide$new(
    utc = tms,
    latitude = 52.3868,
    longitude = 9.7144,
    elevation = 110,
    gravity = 9.8127,
    cutoff = 1.0e-10,
    catalog = "hw95s",
    wave_groups = wave_groups
  )


  tmp1 <- et$predict(method = "gravity", astro_update = 1L)
  tmp300 <- et$predict(method = "gravity", astro_update = 300L)
  expect_equal(tmp1$tide(), tmp300$tide())


  tmp2 <- et$predict(method = "gravity", n_thread = 2)
  expect_equal(tmp1$tide(), tmp2$tide())



})



tms <- as.POSIXct("1990-01-01", tz = "UTC") + 0:86400

wave_groups <- data.frame(start = 0, end = 8)

et <- Earthtide$new(
  utc = tms,
  latitude = 52.3868,
  longitude = 9.7144,
  elevation = 110,
  gravity = 9.8127,
  cutoff = 1.0e-10,
  catalog = "hw95s",
  wave_groups = wave_groups
)

system.time(
tmp <- et$predict(method = "gravity", n_thread = 10, astro_update = 300)
)


tms10 <- as.POSIXct("1990-01-01", tz = "UTC") + seq(0, 86400, 10)

wave_groups <- data.frame(start = 0, end = 8)

et <- Earthtide$new(
  utc = tms10,
  latitude = 52.3868,
  longitude = 9.7144,
  elevation = 110,
  gravity = 9.8127,
  cutoff = 1.0e-10,
  catalog = "hw95s",
  wave_groups = wave_groups
)

system.time({
    et$predict(method = "gravity", n_thread = 10, astro_update = 30)
    tmp <- et$interpolate(tms)
})

tide <- et$tide()




