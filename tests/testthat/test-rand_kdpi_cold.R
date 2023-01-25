test <- rand_kdpi_cold(n=25)

test_that("rand_cold_kdpi() returns right size", {
  expect_length(test[,1], 25)
})

test_that("rand_cold_kdpi() cold iscemia times >= 0", {
  expect_true(all(test[,1])>=0)
})

test_that("rand_cold_kdpi() KDPI >=0 & <=100", {
  expect_true(all(test[,2]>=0 & test[,2]<=100))
})
