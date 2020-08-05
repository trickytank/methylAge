test_that("No errors or warnings", {
  expect_silent(phage <- phenoage(met))
})

phage <- phenoage(met)
# dput(as.data.frame(phage))

phenoage_expected <- structure(list(
  ID = c("GSM1871369", "GSM1871370", "GSM1871371",
         "GSM1871372", "GSM1871373", "GSM1871374", "GSM1871375", "GSM1871376",
         "GSM1871377", "GSM1871378"),
  phenoage = c(18.562589291984, 9.50990252386335,
               -0.829092100213224, -4.97180647309536, -27.7557605016348, 7.43154535791978,
               3.2147952705816, 4.95222803527529, 25.2813249259989, 30.0778049122663)
  ),
  row.names = c(NA, -10L),
  class = "data.frame")

test_that("Expected result of Phenoage", {
  expect_equal(as.data.frame(phage), phenoage_expected)
})
