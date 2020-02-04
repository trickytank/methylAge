check_methylation_data <- function(x, dim_warning = TRUE) {
  assert_true(is.matrix(x) || is.data.frame(x))
  assert_false(is.null(rownames(x)))
  assert_false(is.null(colnames(x)))
  if(dim_warning && nrow(x) > ncol(x)) {
    warning("Number of samples (rows) is more than the number of probes (rows)", immediate. = TRUE)
  }
}
