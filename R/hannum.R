#' Hannum DNA methylation age calculation
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
hannum_clock <- function(x, id_col = "ID", age_col = "hannum_mage", allow_missing = FALSE, dim_warning = TRUE) {
  # Calculate Hannum Methylation Age
  generic_clock(
                x, coef = hannum_coef,
                id_col = id_col, age_col = age_col,
                intercept_name = NULL,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Hannum")
}
