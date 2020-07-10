#' Hannum DNA methylation age calculation
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
hannum_clock <- function(x, id_col = "ID", age_col = "hannum_mAge", allow_missing = FALSE, dim_warning = TRUE) {
  # Calculate Hannum Methylation Age
  check_methylation_data(x, dim_warning = dim_warning)

  # Load coefficients
  coefs_clock <- setNames(hannum_coef$Coefficient, hannum_coef$Marker)

  # Check for missing coefficients
  probes_check <- names(coefs_clock) %in% row.names(x)
  if(!any(probes_check)) {
    stop("No probes from the Hannum clock are in the data.")
  }
  if(allow_missing) {
    probes_present <- names(coefs_clock)[probes_check]
    coefs_clock <- coefs_clock[probes_present]
  } else {
    if(!all(probes_check)) {
      stop("Some Hannum probes are missing from the data.")
    }
    probes_present <- names(coefs_clock)
  }

  # Calculate Hannum Methylation Age
  x_mat <- as.matrix( x[probes_present, ] )
  m_age <- coefs_clock %*% x_mat
  # Cleanup
  tibble::tibble(!!id_col := colnames(m_age), !!age_col := m_age[1,])
}
