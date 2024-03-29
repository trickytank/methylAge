#' Hannum DNA methylation age calculation
#'
#' Calculate the \insertCite{HANNUM2013359;textual}{methylAge} DNA methylation age.
#' Please cite the referenced article if using this function.
#'
#' @inheritParams generic_clock
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
hannum_clock <- function(x, id_out = "ID", age_out = "hannum_mage",
                         allow_missing = getOption('methylAge.allow_missing'),
                         dim_warning = getOption('methylAge.dim_warning')) {
  # Calculate Hannum Methylation Age
  generic_clock(
                x, coef = hannum_coef,
                id_out = id_out, age_out = age_out,
                intercept_name = NULL,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Hannum")
}
