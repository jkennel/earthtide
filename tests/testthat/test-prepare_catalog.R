context("test-prepare_catalog")

test_that("prepare_catalog works", {
  
  expect_warning(.prepare_catalog(1e-10, NA, catalog = 'ksm04'))
    
})
