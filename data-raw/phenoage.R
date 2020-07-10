# ---- Prepare PhenoAge Coefficients ----
# Supplementary files will need to be obtained from the publication listed below, and
# placed in the data-raw/supplementary/ directory

# # PhenoAge clock:
# Requires supplement 2, Table S6 (aging-10-101414-s002.csv) from:
# Levine ME, Lu AT, Quach A, et al.
# An epigenetic biomarker of aging for lifespan and healthspan.
# Aging (Albany NY). 2018;10(4):573-591.
# doi:10.18632/aging.101414
# https://www.aging-us.com/article/101414/text

library(devtools)
library(tidyverse)
library(tools)

local_phenoage <- "data-raw/supplementary/PhenoAge/aging-10-101414-s002.csv"

if(md5sum(local_phenoage) != "8a812be441c2d1c553bb6cd97597562d") {
  stop("PhenoAge coefficients file does not appear to be correct according to the md5sum.")
}

phenoage_coef_raw <- readr::read_csv(local_phenoage)
# Make a syntactically valid name
phenoage_coef <- phenoage_coef_raw %>% rename_all(~ str_replace_all(.x, " ", "_"))

use_data(phenoage_coef)
