#' Phenoage DNA methylation age calculation
#'
#' Calculate the \insertCite{phenoage;textual}{methylAge} DNA methylation Phenoage.
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
#' @importFrom Rdpack reprompt
phenoage <- function(x, id_col = "ID", age_col = "phenoage",
                     allow_missing = getOption('methylAge.allow_missing'),
                     dim_warning = getOption('methylAge.dim_warning')) {
  # Calculate Phenoage Methylation Age
  generic_clock(
                x, coef = phenoage_coef,
                id_col = id_col, age_col = age_col,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Phenoage"
  )
}
