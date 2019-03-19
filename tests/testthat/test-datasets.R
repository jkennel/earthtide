context("test-datasets")

# test_that("get dut1 works", {
# 
#   obs <- get_dut1()
#   
#   expect_s3_class(obs, 'data.frame')
#   expect_equal(ncol(obs), 9L)
#   expect_gt(nrow(obs), 20000)
#   
# })


test_that("get iers works", {
  
  dut1 <- get_dut1_iers()
  
  expect_s3_class(dut1, 'data.frame')
  expect_equal(ncol(dut1), 8L)
  expect_gt(nrow(dut1), 20000)
  
})
