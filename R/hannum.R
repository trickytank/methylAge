#' @export
#' @import checkmate
hannum <- function(x, id_col = "ID", age_col = "hannum_mAge", dim_warning = TRUE) {
  # Calculate Hannum Methylation Age
  check_methylation_data(x, dim_warning = dim_warning)

  # Load coefficients
  coefs_clock <- setNames(hannum_coef$Coefficient, hannum_coef$Marker)
  # Calculate Hannum Methylation Age
  x_mat <- as.matrix( x[, names(coefs_clock)] )
  m_age <- x_mat %*% coefs_clock
  # Cleanup
  colnames(m_age) <- age_col
  m_age_tbl <- as_tibble(m_age, rownames = id_col)
  m_age_tbl
}
