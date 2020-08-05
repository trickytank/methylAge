test_that("No errors or warnings", {
  expect_silent(horvath_age <- horvath_clock(met))
})

horvath_age <- horvath_clock(met)
# dput(as.data.frame(horvath_age))

horvath_expected <- structure(list(
  ID = c("GSM1871369", "GSM1871370", "GSM1871371",
         "GSM1871372", "GSM1871373", "GSM1871374", "GSM1871375", "GSM1871376",
         "GSM1871377", "GSM1871378"),
  horvath_mage = c(14.1429002501734,
                   8.10135386746589, 5.01534231478755, 4.11179252376334, 2.04644977773932,
                   7.37981277129947, 4.98860349899741, 6.20351880191709, 25.4473268085755,
                   34.8953324503867)
  ),
  row.names = c(NA, -10L),
  class = "data.frame")

test_that("Expected result of horvath clock", {
  expect_equal(as.data.frame(horvath_age), horvath_expected)
})
