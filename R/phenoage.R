#' Hannum DNA methylation age calculation
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
phenoage_clock <- function(x, id_col = "ID", age_col = "phenoage", allow_missing = FALSE, dim_warning = TRUE) {
  # Calculate Hannum Methylation Age
  generic_clock(
                x, coef = phenoage_coef,
                marker_col = "CpG", coef_col = "Weight",
                id_col = id_col, age_col = age_col,
                allow_missing = allow_missing, dim_warning = dim_warning,
                clock_name = "Phenoage"
  )
}
