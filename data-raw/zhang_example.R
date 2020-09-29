# --- Download Zhang example data
# file is downloaded automatically.
#
# Zhang Q, Vallerga C, Walker R, Lin T, Henders A, Montgomery G, He J, Fan D, Fowdar J, Kennedy M, et al:
# Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing.
# Genome Medicine, 2019.
#
# https://github.com/qzhang314/DNAm-based-age-predictor

library(devtools)
library(tidyverse)
library(tools)
library(fs)

local_zhang_methylation <- 'data-raw/supplementary/zhang/data.rds'
local_zhang_age <- 'data-raw/supplementary/zhang/data.age'

dir_create(path_dir(local_zhang_methylation))

if(!file.exists(local_zhang_methylation)) {
  download.file('https://raw.githubusercontent.com/qzhang314/DNAm-based-age-predictor/13b213b9644d6661556cf30129f1c17a22c57274/data.rds', local_zhang_methylation)
}
if(md5sum(local_zhang_methylation) != "b5af935a6fbb89a6dc632243e45c69da") {
  stop("Zhang methylation example data file does not appear to be correct according to the md5sum.")
}

if(!file.exists(local_zhang_age)) {
  download.file('https://raw.githubusercontent.com/qzhang314/DNAm-based-age-predictor/13b213b9644d6661556cf30129f1c17a22c57274/data.age', local_zhang_age)
}
if(md5sum(local_zhang_age) != "81965c06e7b85e8c8f98928244e34b5e") {
  stop("Zhang methylation example age file does not appear to be correct according to the md5sum.")
}

met <- t(readRDS(local_zhang_methylation))

met_age <- readr::read_delim(local_zhang_age, " ")

use_data(met)
use_data(met_age)
