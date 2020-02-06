#' Hannum DNA methylation age calculation
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
hannum <- function(x, id_col = "ID", age_col = "hannum_mAge", dim_warning = TRUE) {
  # Calculate Hannum Methylation Age
  check_methylation_data(x, dim_warning = dim_warning)

  # Load coefficients
  coefs_clock <- setNames(hannum_coef$Coefficient, hannum_coef$Marker)
  # Calculate Hannum Methylation Age
  x_mat <- as.matrix( x[names(coefs_clock), ] )
  m_age <- coefs_clock %*% x_mat
  # Cleanup
  tibble::tibble(!!id_col := colnames(m_age), !!age_col := m_age[1,])
}
