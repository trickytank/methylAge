#' Lu DNA methylation telomere length (DNAmTL) calculation
#'
#' @inheritParams generic_clock
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
lu_telomere <- function(x, id_out = "ID", age_out = "DNAmTL",
                        allow_missing = getOption('methylAge.allow_missing'),
                        dim_warning = getOption('methylAge.dim_warning')) {
  # Calculate Hannum Methylation Age
  generic_clock(
                x, coef = lu_coef,
                id_out = id_out, age_out = age_out,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Lu DNAmTL"
  )
}
