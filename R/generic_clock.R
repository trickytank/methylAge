#' Generic linear model DNA methylation age calculation
#'
#' This function gives a methylation age based on the coefficients of a trained clock that
#' does not transform the data. Transforms may be performed outside of this function.
#' Several clocks in this package use this function for methylation age calculation.
#'
#' @param x Matrix or data.frame of methylation proportions, with probes/markers as `rownames` and samples as `colnames`.
#' @param coef data.frame of coefficients, with a column for marker and column for model coefficient.
#' @param col_marker Name of the marker/probe column in `coef`.
#' @param col_coef Name of the coefficient column in `coef`
#' @param intercept_name If character, the name of the intercept marker/probe in `coef`.
#'                       If numeric, then this is the intercept value.
#'                       If NULL, then the intercept is set to 0.
#' @param id_out Name of the sample ID in the output tibble.
#' @param age_out Name of the estimated age column in the output tibble.
#' @param allow_missing If FALSE, then all markers in `coef` must be present in `x.`
#' @param dim_warning If TRUE, then a warning is raised when there are more columns than rows in `x`.
#' @param clock_name Name of clock for errors and warnings.
#'
#' @return A tibble with and id column and estimated methylation age.
#'
#' @examples
#'
#' # Run the Hannum clock.
#' generic_clock(met, coef = as.data.frame(hannum_coef), intercept_name = NULL, clock_name = "Hannum")
#'
#' # Run the Phenoage clock.
#' generic_clock(met, coef = phenoage_coef, clock_name = "Phenoage")
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
#' @importFrom stats setNames
generic_clock <- function(x, coef,
                          col_marker = "marker", col_coef = "coefficient",
                          intercept_name = "Intercept",
                          id_out = "ID", age_out = "mage",
                         allow_missing = getOption('methylAge.allow_missing'),
                         dim_warning = getOption('methylAge.dim_warning'),
                         clock_name = "generic"
                         ) {
  # Calculate generic linear model DNA methylation age
  check_methylation_data(x, dim_warning = dim_warning)

  # Load coefficients
  coefs_clock <- setNames(coef[[col_coef]], coef[[col_marker]])

  # Extract intercept if present
  if(!is.null(intercept_name)) {
    if(is.numeric(intercept_name)) {
      intercept <- intercept_name
    } else if(is.character(intercept_name)) {
      intercept <- coefs_clock[intercept_name]
      coefs_clock <- coefs_clock[names(coefs_clock) != intercept_name]
      if(is.na(intercept)) {
        stop("Value for intercept with name '", intercept_name, "' not found in coef.")
      }
      assert_number(intercept, finite = TRUE)
    } else {
      stop("intercept name is not character or numeric.")
    }
  } else {
    intercept <-  0
  }

  # Check for missing coefficients
  probes_check <- names(coefs_clock) %in% row.names(x)
  if(!any(probes_check)) {
    stop("No probes from the ", clock_name, " estimator is in the data.")
  }
  if(allow_missing) {
    probes_present <- names(coefs_clock)[probes_check]
    coefs_clock <- coefs_clock[probes_present]
  } else {
    if(!all(probes_check)) {

      stop(sum(!probes_check), " ", clock_name, " probes are missing from the data.")
    }
    probes_present <- names(coefs_clock)
  }

  # Calculate Methylation Age
  x_mat <- as.matrix( x[probes_present, ] )
  m_age <- coefs_clock %*% x_mat + intercept
  # Cleanup
  tibble::tibble(!!id_out := colnames(m_age), !!age_out := m_age[1,])
}
