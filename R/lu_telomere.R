#' Lu DNA methylation telomere length (DNAmTL) calculation
#'
#' Calculate the \insertCite{lu_telomere;textual}{methylAge} DNA methylation age.
#' Please cite the referenced article if using this function.
#'
#' @inheritParams generic_clock
#' @param est_out Name of the estimated telomere length column in the output tibble.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
lu_telomere <- function(x, id_out = "ID", est_out = "DNAmTL",
                        allow_missing = getOption('methylAge.allow_missing'),
                        dim_warning = getOption('methylAge.dim_warning')
                        ) {
  # Calculate Hannum Methylation Age
  generic_clock(
                x, coef = lu_coef,
                id_out = id_out, age_out = est_out,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Lu DNAmTL"
  )
}
