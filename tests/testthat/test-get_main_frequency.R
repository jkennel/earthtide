test_that("get_main_frequency works", {
  
  
  expect_equal(get_main_frequency(c(0.89,1.8), c(0.94,2)),
               c(0.9295357, 1.9322736))
  
  expect_equal(get_main_frequency(0, 6),
               1.9322736)
  
  expect_error(get_main_frequency(0, c(2,1)))
  expect_error(get_main_frequency(2, 1))
  
  
})
