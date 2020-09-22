#' Convenience function for DunedinPoAm pace of ageing
#'
#' Calculate the (TODO: cite https://doi.org/10.7554/eLife.54870) cortex tissueDNA methylation age calculation.
#' Please cite the referenced article if using this function.
#'
#' This relies on installing the package with `
#' ```
#' devtools::install_github("danbelsky/DunedinPoAm38")
#' ```
#'
#' @inheritParams generic_clock
#' @param proportionOfProbesRequired The threshold for missing data, as specified in DunedinPoAm38::PoAmProjector.
#'
#' @references
#' \insertAllCited{}
#'
#' @export
#' @import checkmate
#' @import dplyr
#' @import tibble
#' @importFrom Rdpack reprompt
dunedinpoam_helper <- function(x,
                     proportionOfProbesRequired = 0.8,
                     id_out = "ID", est_out = "DunedinPoAm",
                     allow_missing = getOption('methylAge.allow_missing'),
                     dim_warning = getOption('methylAge.dim_warning')) {
  # Calculate DunedinPoAm
  if("DunedinPoAm38" %in% rownames(installed.packages())) {
    check_methylation_data(x, dim_warning = dim_warning)

    # TODO: check allow missing
    probes_check <- DunedinPoAm38::mPOA_Models$model_probes$DunedinPoAm_38 %in% row.names(x)
    if(!any(probes_check)) {
      stop("No probes from the DunedinPoAm estimator is in the data.")
    }
    if(!allow_missing && !all(probes_check)) {
      stop(sum(!probes_check), " DunedinPoAm probes are missing from the data.")
    }

    poAm <- DunedinPoAm38::PoAmProjector(x, proportionOfProbesRequired)
    dunedinpoam_clean(poAm, id_out, est_out)

  } else {
    stop("The package DunedinPoAm38 should be installed from https://github.com/danbelsky/DunedinPoAm38 to run this function.")
  }
}

#' Clean DunedinPoAm output to fit with the output of other estimators in the MethylAge package
dunedinpoam_clean <- function(x, id_out = "ID", est_out = "DunedinPoAm") {
  y <- tibble(
    id_pre = names(x[[1]]),
    est_out_pre = x[[1]]
  )
  names(y) <- c(id_out, est_out)
  y
}
