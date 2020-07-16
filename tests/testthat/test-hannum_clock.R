test_that("No errors or warnings", {
  expect_silent(hannum_age <- hannum_clock(met))
})

hannum_age <- hannum_clock(met)

hannum_expected <- structure(list(
  ID = c("GSM1871369", "GSM1871370", "GSM1871371",
         "GSM1871372", "GSM1871373", "GSM1871374", "GSM1871375", "GSM1871376",
         "GSM1871377", "GSM1871378"),
  hannum_mage = c(27.260503326069,
                  18.6222871384364, 12.8078073218327, 9.04056419859078, 0.819221894980454,
                  17.8833320203084, 12.9488556518456, 15.0877531138839, 36.4150319915383,
                  43.6277412999478)),
  row.names = c(NA, -10L),
  class = "data.frame")

test_that("Expected result of hannum clock", {
  expect_equal(as.data.frame(hannum_age), hannum_expected)
})
