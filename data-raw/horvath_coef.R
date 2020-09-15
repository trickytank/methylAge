# --- Prepare Horvath Coefficients data
# Coefficient file is downloaded automatically.
# If this does not work automatically, then download 'AdditionalFile3.csv' at
# https://horvath.genetics.ucla.edu/html/dnamage/
# and store in the location specified by local_horvath below.
# In the original publication, additional file 3 has a description line at the top;
# if using this file, delete the first two lines. For example with the Linux command:
# tail -n +3 gb-2013-14-10-r115-S3.csv > AdditionalFile3.csv
# (supplemetary file name may vary depending on source)

# Horvath S. DNA methylation age of human tissues and cell types. Genome Biol.
# 2013;14(10):R115.
# PubMed PMID: 24138928; PubMed Central PMCID: PMC4015143.
# https://doi.org/10.1186/gb-2013-14-10-r115

library(devtools)
library(tidyverse)
library(tools)
library(fs)

local_horvath <- 'data-raw/supplementary/horvath/AdditionalFile3.csv'

dir_create(path_dir(local_horvath))

if(!file.exists(local_horvath)) {
  download.file('https://horvath.genetics.ucla.edu/html/dnamage/AdditionalFile3.csv', local_horvath)
}
if(md5sum(local_horvath) != "27df5691e91ee6dc362458efaae1a739") {
  stop("Horvath coefficients file does not appear to be correct according to the md5sum.")
}

horvath_coef_raw <- read_csv(local_horvath)
horvath_coef <- horvath_coef_raw %>%
  select(marker = CpGmarker, coefficient = CoefficientTraining) %>%
  mutate(marker = if_else(marker == "(Intercept)", "Intercept", marker))

use_data(horvath_coef)


# And Horvath annotations
local_probe_21kdat <- 'data-raw/supplementary/horvath/probeAnnotation21kdatMethUsed.csv'
local_probe_27k <- 'data-raw/supplementary/horvath/datMiniAnnotation27k.csv'

if(!file.exists(local_probe_21kdat)) {
  download.file('https://horvath.genetics.ucla.edu/html/dnamage/probeAnnotation21kdatMethUsed.csv', local_probe_21kdat)
}
if(md5sum(local_probe_21kdat) != "f8eb4f868e328f4c4f492bbf83ea0c4d") {
  stop("Horvath annotation 21k datMethUsed file does not appear to be correct according to the md5sum.")
}

if(!file.exists(local_probe_27k)) {
  download.file('https://horvath.genetics.ucla.edu/html/dnamage/datMiniAnnotation27k.csv', local_probe_27k)
}
if(md5sum(local_probe_27k) != "bad8eed11b7d3ad520a4b141cc393f29") {
  stop("Horvath mini annotation 27k file does not appear to be correct according to the md5sum.")
}

horvath_probeAnnotation21kdatMethUsed <- read.csv(local_probe_21kdat)
horvath_probeAnnotation27k <- read.csv(local_probe_27k)

use_data(horvath_probeAnnotation21kdatMethUsed)
use_data(horvath_probeAnnotation27k)
