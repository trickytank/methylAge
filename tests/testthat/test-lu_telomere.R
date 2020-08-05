test_that("No errors or warnings", {
  expect_silent(telomere_est <- lu_telomere(met))
})

telomere_est <- lu_telomere(met)
# dput(as.data.frame(telomere_est))

telomere_expected <- structure(list(
  ID = c("GSM1871369", "GSM1871370", "GSM1871371",
         "GSM1871372", "GSM1871373", "GSM1871374", "GSM1871375", "GSM1871376",
         "GSM1871377", "GSM1871378"),
  DNAmTL = c(7.66646738654308, 7.93317104396237,
             7.81795202035524, 8.15263996282507, 8.58446824258146, 8.01109994162254,
             8.09002859037805, 8.15199173183434, 7.52191563222574, 6.99250190349181)
  ),
  row.names = c(NA, -10L),
  class = "data.frame"
)

test_that("Expected result of DNA methylation Telomere estimate", {
  expect_equal(as.data.frame(telomere_est), telomere_expected)
})
