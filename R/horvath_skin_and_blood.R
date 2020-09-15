#' Phenoage DNA methylation age calculation
#'
#' Calculate the \insertCite{horvath_skin_blood;textual}{methylAge} Skin and Blood methylation age.
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
horvath_skin_blood_clock <- function(x, id_out = "ID", age_out = "horvath_skin_blood",
                     allow_missing = getOption('methylAge.allow_missing'),
                     dim_warning = getOption('methylAge.dim_warning')) {
  # Calculate Horvath Skin and Blood Methylation Age
  mage <- generic_clock(
    x, coef = hsb_coef,
    id_out = id_out, age_out = age_out,
    allow_missing = allow_missing, dim_warning = dim_warning,
    clock_name = "Horvath skin and blood"
  )
  mage[[2]] <- anti.trafo(mage[[2]])
  mage
}
