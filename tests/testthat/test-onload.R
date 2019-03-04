context("zzz.R - .onLoad")

test_that(".onLoad loads data", {
  
  obj_names <- c("hw95s", "dut1", "simon_coef_1994", 'ksm04')

  expect_true(all(obj_names %in% data(package ='earthtide')$results[,3]))
  
})
