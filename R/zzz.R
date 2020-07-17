# Package functions and documentation

methylAge_default_options <- list(
  methylAge.allow_missing = FALSE,
  methylAge.dim_warning = TRUE
)

.onLoad <- function(libname, pkgname) {
  op <- options()
  toset <- !(names(methylAge_default_options) %in% names(op))
  if (any(toset)) options(methylAge_default_options[toset])

  invisible()
}
