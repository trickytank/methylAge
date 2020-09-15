# ---- Prepare Horvath's DNAm Skin and Blood clock Coefficients ----
# Supplementary files will need to be obtained from the publication listed below, and
# placed in the data-raw/supplementary/ directory

# # Lu DNAmTL estimator
# Requires Supplementary Dataset 2 (JHuEc465v95mjTR55_sd5.csv) from:
# Horvath S, Oshima J, Martin GM, Lu AT, Quach A, Cohen H, Felton S, Matsuyama M, Lowe D, Kabacik S, Wilson JG, Reiner AP, Maierhofer A, et al, .
# Epigenetic clock for skin and blood cells applied to Hutchinson Gilford Progeria Syndrome and ex vivo studies.
# Aging (Albany NY). 2018; 10:1758-1775.
# https://doi.org/10.18632/aging.101508

library(devtools)
library(tidyverse)
library(tools)

local_hsb <- "data-raw/supplementary/horvath_skin_and_blood/JHuEc465v95mjTR55_sd5.csv"

if(tools::md5sum(local_hsb) != "a06862d6b6ff91ad969563d3f44716c2") {
  stop("Horvath Skin and Blood DNA methylation clock coefficients file does not appear to be correct according to the md5sum.")
}

hsb_coef_raw <- readr::read_csv(local_hsb)
# Make a syntactically valid name
hsb_coef <- hsb_coef_raw %>%
  select(marker = ID, coefficient = Coef) %>%
  mutate(marker = if_else(marker == "(Intercept)", "Intercept", marker))

use_data(hsb_coef)
