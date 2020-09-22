#' Shireby's cortex tissue DNA methylation age calculation
#'
#' Calculate the \insertCite{Shireby2020_cortex;textual}{methylAge} cortex tissueDNA methylation age calculation.
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
shireby_cortex_clock <- function(x, id_out = "ID", age_out = "shireby_cortex_mage",
                     allow_missing = getOption('methylAge.allow_missing'),
                     dim_warning = getOption('methylAge.dim_warning')) {
  # Calculate Shireby cortex tissue clock.
  mage <- generic_clock(
    x, coef = shireby_coef,
    id_out = id_out, age_out = age_out,
    allow_missing = allow_missing, dim_warning = dim_warning,
    clock_name = "Shireby cortex tissue DNA methylation age clock"
  )
  mage[[2]] <- anti.trafo(mage[[2]])
  mage
}
