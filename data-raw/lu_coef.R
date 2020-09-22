# ---- Prepare Lu's DNAmTL Coefficients ----
# Supplementary files will need to be obtained from the publication listed below, and
# placed in the data-raw/supplementary/ directory

# # Lu DNAmTL estimator
# Requires Supplementary Data 1 (aging-11-102173-s003.xlsx) from:
# Lu AT, Seeboth A, Tsai PC, et al.
# DNA methylation-based estimator of telomere length.
# Aging (Albany NY). 2019;11(16):5895-5923.
# doi:10.18632/aging.102173
# https://dx.doi.org/10.18632%2Faging.102173

library(devtools)
library(tidyverse)
library(tools)
library(readxl)

local_lu <- "data-raw/supplementary/lu_DNAmTL/aging-11-102173-s003.xlsx"

if(md5sum(local_lu) != "ccf1be3cd87ce10f1c5b58213b0fafdd") {
  stop("Lu DNAmTL coefficients file does not appear to be correct according to the md5sum.")
}

lu_coef_raw <- readxl::read_excel(local_lu, skip = 5)
# Make consistent across package
lu_coef <- lu_coef_raw %>%
  select(marker = Variable, coefficient = Coefficient)

use_data(lu_coef)
