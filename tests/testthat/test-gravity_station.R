context("test-gravity_station")

test_that("gravity_station works", {

  expect_equal(gravity_station(0,0), 9.78032677, tolerance = 1e-6)
  expect_equal(gravity_station(90, 0), 9.832186, tolerance = 1e-6)
  expect_equal(gravity_station(45, 0), 9.806199, tolerance = 1e-6)

  expect_equal(gravity_station(0,1000), 9.777257, tolerance = 1e-6)
  expect_equal(gravity_station(0, 10000), 9.749696, tolerance = 1e-6)

})
