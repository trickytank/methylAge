# Smaller dataset to test with :p
metm <- met[1:4, 1:3]

test_that("No warnings or errors on clean data", {
  expect_silent(check_methylation_data(metm))
})

test_that("Warning when more rows than columns", {
  expect_warning(check_methylation_data(metm[1:2, 1:3]))
})

test_that("Error with wrong data type", {
  expect_error(check_methylation_data(1))
})

metm_col <- metm
metm_row <- metm
colnames(metm_col) <- rownames(metm_row) <- NULL
test_that("Error when no column or row names", {
  expect_error(check_methylation_data(metm_col))
  expect_error(check_methylation_data(metm_row))
})
