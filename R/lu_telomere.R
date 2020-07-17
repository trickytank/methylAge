#' Lu DNA methylation telomere length (DNAmTL) calculation
#'
#' @inheritParams generic_clock
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
lu_telomere <- function(x, id_col = "ID", age_col = "DNAmTL", allow_missing = FALSE, dim_warning = TRUE) {
  # Calculate Hannum Methylation Age
  generic_clock(
                x, coef = lu_coef,
                id_col = id_col, age_col = age_col,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Lu DNAmTL"
  )
}
