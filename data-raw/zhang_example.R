# --- Download Zhang example data
# file is downloaded automatically.
#

#TODO work in progress

# Zhang Q, Vallerga C, Walker R, Lin T, Henders A, Montgomery G, He J, Fan D, Fowdar J, Kennedy M, et al:
# Improved precision of epigenetic clock estimates across tissues and its implication for biological ageing.
# Genome Medicine, 2019.

library(devtools)
library(tidyverse)
library(tools)
library(fs)

local_zhang_methylation <- 'data-raw/supplementary/zhang/data.rds'


#
dir_create(path_dir(local_zhang_methylation))

if(!file.exists(local_zhang_methylation)) {
  download.file('https://raw.githubusercontent.com/qzhang314/DNAm-based-age-predictor/13b213b9644d6661556cf30129f1c17a22c57274/data.rds', local_zhang_methylation)
}
if(md5sum(local_zhang_methylation) != "b5af935a6fbb89a6dc632243e45c69da") {
  stop("Zhang methylation example data file does not appear to be correct according to the md5sum.")
}

met <- t(readRDS(local_zhang_methylation))

use_data(met)
