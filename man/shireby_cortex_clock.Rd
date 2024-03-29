% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/shireby_cortex.R
\name{shireby_cortex_clock}
\alias{shireby_cortex_clock}
\title{Shireby's cortex tissue DNA methylation age calculation}
\usage{
shireby_cortex_clock(
  x,
  id_out = "ID",
  age_out = "shireby_cortex_mage",
  allow_missing = getOption("methylAge.allow_missing"),
  dim_warning = getOption("methylAge.dim_warning")
)
}
\arguments{
\item{x}{Matrix or data.frame of methylation proportions, with probes/markers as `rownames` and samples as `colnames`.}

\item{id_out}{Name of the sample ID in the output tibble.}

\item{age_out}{Name of the estimated age column in the output tibble.}

\item{allow_missing}{If FALSE, then all markers in `coef` must be present in `x.`}

\item{dim_warning}{If TRUE, then a warning is raised when there are more columns than rows in `x`.}
}
\description{
Calculate the \insertCite{Shireby2020_cortex;textual}{methylAge} cortex tissueDNA methylation age calculation.
Please cite the referenced article if using this function.
}
\references{
\insertAllCited{}
}
