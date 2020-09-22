# ---- Prepare Shireby's DNAmAge Coefficients for human cortex tissue ----
# Supplementary files will need to be obtained from the publication listed below, and
# placed in the data-raw/supplementary/ directory

# # Shireby's DNAmAge estimator for human cortex tissue
# Requires Supplementary Table 2 in 063719_file06.xlsx from:
# Gemma L Shireby, Jonathan P Davies, Paul T Francis, Joe Burrage, Emma M Walker, Grant W A Neilson, Aisha Dahir, Alan J Thomas, Seth Love, Rebecca G Smith, Katie Lunnon, Meena Kumari, Leonard C Schalkwyk, Kevin Morgan, Keeley Brookes, Eilis J Hannon, Jonathan Mill
# *Recalibrating the Epigenetic Clock: Implications for Assessing Biological Age in the Human Cortex*
# bioRxiv 2020.04.27.063719; doi: https://doi.org/10.1101/2020.04.27.063719

library(devtools)
library(tidyverse)
library(tools)
library(readxl)

local_shireby <- "data-raw/supplementary/shireby/063719_file06.xlsx"

if(md5sum(local_shireby) != "d9434f65b2855131c1a8f8b117a25fcd") {
  stop("Shireby DNAmAge coefficients file does not appear to be correct according to the md5sum.")
}

shireby_coef_raw <- readxl::read_excel(local_shireby, sheet = "Supplementary Table 2")
# Make consistent across package
shireby_coef <- shireby_coef_raw %>%
  select(marker = probe, coefficient = coef) %>%
  mutate(marker = if_else(marker == "intercept", "Intercept", marker))

use_data(shireby_coef)
