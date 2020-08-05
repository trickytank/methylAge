test_that("No errors or warnings", {
  expect_message(expect_warning(zhang_age <- zhang_clock(met), regexp = NA))
})

zhang_age <- zhang_clock(met)
# dput(as.data.frame(zhang_age))

zhang_expected <- structure(list(
  ID = c("GSM1871369", "GSM1871370", "GSM1871371",
         "GSM1871372", "GSM1871373", "GSM1871374", "GSM1871375", "GSM1871376",
         "GSM1871377", "GSM1871378"),
  zhang_en_mage = c(12.8741932290722,
                    8.6265478302361, 1.45133392533948, -3.74265502664618, -7.6513650649452,
                    6.44173638841287, 0.267038063468874, 3.4376646456282, 26.6353286892498,
                    38.0091473731159),
  zhang_blup_mage = c(13.524130747235, 9.86968238281824,
                      6.6422769698528, 1.99624082557108, 1.73798593075117, 8.51620943636564,
                      6.07119211359762, 7.95060847098556, 23.4502284263066, 34.8608125721457)
  ),
  row.names = c(NA, -10L),
  class = "data.frame")

test_that("Expected result of zhang clock", {
  expect_equal(as.data.frame(zhang_age), zhang_expected)
})
