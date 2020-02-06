#' Check methylation data is formatted in a way suitable for the methylAge package
#' @export
#' @import checkmate
check_methylation_data <- function(x, dim_warning = TRUE) {
  assert_true(is.matrix(x) || is.data.frame(x))
  assert_false(tibble::is_tibble(x))
  assert_false(is.null(rownames(x)))
  assert_false(is.null(colnames(x)))
  if(dim_warning && nrow(x) < ncol(x)) {
    warning("Number of samples (columns) is more than the number of probes (rows)", immediate. = TRUE)
  }
}
