# --- Prepare Zhang Coefficients data
# Coefficient file is downloaded automatically.
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

local_zhang_en <- 'data-raw/supplementary/zhang/en.coef'
local_zhang_blup <- 'data-raw/supplementary/zhang/blup.coef'

dir_create(path_dir(local_zhang_en))

if(!file.exists(local_zhang_en)) {
  download.file('https://raw.githubusercontent.com/qzhang314/DNAm-based-age-predictor/master/en.coef', local_zhang_en)
}
if(md5sum(local_zhang_en) != "4852d2fc5ddd348eb94bca279109c546") {
  stop("Zhang Elastic Net coefficients file does not appear to be correct according to the md5sum.")
}

if(!file.exists(local_zhang_blup)) {
  download.file('https://raw.githubusercontent.com/qzhang314/DNAm-based-age-predictor/master/blup.coef', local_zhang_blup)
}
if(md5sum(local_zhang_blup) != "11a8e62991b1b77566bb10e25da0a48c") {
  stop("Zhang BLUP coefficients file does not appear to be correct according to the md5sum.")
}

zhang_en_coef <- read.table(local_zhang_en, stringsAsFactor = FALSE, header = TRUE)
zhang_blup_coef <- read.table(local_zhang_blup, stringsAsFactor = FALSE, header = TRUE)

zhang_en_coef <- zhang_en_coef %>% rename(marker = probe, coefficient = coef)
zhang_blup_coef <- zhang_blup_coef %>% rename(marker = probe, coefficient = coef)

use_data(zhang_en_coef, zhang_blup_coef)
