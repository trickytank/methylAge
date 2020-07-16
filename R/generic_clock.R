#' Generic linear model DNA methylation age calculation
#'
#' This function gives a methylation age based on the coefficients of a trained clock that
#' does not transform the data. Transforms may be performed outside of this function.
#' Several clocks in this package use this function for the clock calculation.
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
generic_clock <- function(x, coef,
                          marker_col, coef_col,
                          id_col = "ID", age_col = "mage",
                         allow_missing = FALSE, dim_warning = TRUE,
                         clock_name = "generic") {
  # Calculate generic linear model DNA methylation age
  check_methylation_data(x, dim_warning = dim_warning)

  # Load coefficients
  coefs_clock <- setNames(coef[[coef_col]], coef[[marker_col]])

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
      stop("Some ", clock_name, " probes are missing from the data.")
    }
    probes_present <- names(coefs_clock)
  }

  # Calculate Methylation Age
  x_mat <- as.matrix( x[probes_present, ] )
  m_age <- coefs_clock %*% x_mat
  # Cleanup
  tibble::tibble(!!id_col := colnames(m_age), !!age_col := m_age[1,])
}
