#' Run all methylation age calculations
#'
#' Run all age clocks (or a selection of them).
#' Please cite the referenced articles for each clock if using this function.
#'
#' @inheritParams generic_clock
#'
#' @details
#' This should be used on normalized data.
#'
#' The clocks were developed in
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
all_clocks <- function(x, coef,
                        id_out = "ID",
                        clocks = "age",
                        cols_out = NULL, # age columns out
                        allow_missing = getOption('methylAge.allow_missing'),
                        dim_warning = getOption('methylAge.dim_warning')
) {
  stop("Function not yet written.")
}
