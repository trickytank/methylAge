#' Phenoage DNA methylation age calculation
#'
#' @inheritParams generic_clock
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
phenoage <- function(x, id_col = "ID", age_col = "phenoage",
                     allow_missing = getOption('methylAge.allow_missing'),
                     dim_warning = getOption('methylAge.dim_warning')) {
  # Calculate Hannum Methylation Age
  generic_clock(
                x, coef = phenoage_coef,
                id_col = id_col, age_col = age_col,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Phenoage"
  )
}
