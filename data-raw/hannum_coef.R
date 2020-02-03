# --- Prepare Hannum Coefficients
# Supplementary files will need to be obtained from the publication listed below, and
# placed in the data-raw/supplementary/ directory

# # Hannum clock:
# Requires supplementary Table S3 (1-s2.0-S1097276512008933-mmc2.xlsx) from:
# Gregory Hannum, Justin Guinney, Ling Zhao, Li Zhang, Guy Hughes, SriniVas Sadda, Brandy Klotzle, Marina Bibikova, Jian-Bing Fan, Yuan Gao, Rob Deconde, Menzies Chen, Indika Rajapakse, Stephen Friend, Trey Ideker, Kang Zhang,
# Genome-wide Methylation Profiles Reveal Quantitative Views of Human Aging Rates,
# Molecular Cell,
# Volume 49, Issue 2,
# 2013,
# Pages 359-367,
# ISSN 1097-2765,
# https://doi.org/10.1016/j.molcel.2012.10.016.
# (http://www.sciencedirect.com/science/article/pii/S1097276512008933)

library(devtools)
library(tidyverse)
library(readxl)
library(tools)

local_hannum <- "data-raw/supplementary/1-s2.0-S1097276512008933-mmc2.xlsx"

if(md5sum(local_hannum) != "5b9aeb1ba2ec9915b9f815d9bbc30f9e") {
  stop("Hannum coefficients file does not appear to be correct according to the md5sum.")
}

hannum_coef_raw <- readxl::read_excel(local_hannum, sheet = "Model_PrimaryData")
# Make a syntactically valid name
hannum_coef <- hannum_coef_raw %>% rename(CpG_Island = "CpG Island")

use_data(hannum_coef)
