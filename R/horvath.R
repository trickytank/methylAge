#' Horvath DNA methylation age calculation
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
horvath <- function(x, id_col = "ID", age_col = "horvath_mAge", normalize = FALSE, dim_warning = TRUE) {
  # Calculate Horvath Methylation Age
  check_methylation_data(x, dim_warning = dim_warning)

  if(normalize) {
    x <- bmiq_normalization(x)
  }

  # Load coefficients
  coefs_clock_all <- setNames(horvath_coef$CoefficientTraining, horvath_coef$CpGmarker)
  coefs_clock <- coefs_clock_all[-1] # remove intercept
  # Calculate Horvath Methylation Age
  x_mat <- as.matrix( x[names(coefs_clock), ] )
  m_age <- anti.trafo(coefs_clock %*% x_mat + coefs_clock_all[1])
  # Cleanup
  tibble(!!id_col := colnames(m_age), !!age_col := m_age[1,])
}

trafo <- function(x,adult.age=20) { x=(x+1)/(1+adult.age); y=ifelse(x<=1, log( x),x-1);y }
anti.trafo <- function(x,adult.age=20) { ifelse(x<0, (1+adult.age)*exp(x)-1, (1+adult.age)*x+adult.age) }
