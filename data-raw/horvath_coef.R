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

local_horvath <- 'data-raw/supplementary/AdditionalFile3.csv'

if(!file.exists(local_horvath)) {
  download.file('https://horvath.genetics.ucla.edu/html/dnamage/AdditionalFile3.csv', local_horvath)
}
if(md5sum(local_horvath) != "27df5691e91ee6dc362458efaae1a739") {
  stop("Horvath coefficients file does not appear to be correct according to the md5sum.")
}

horvath_coef_raw <- read.csv(local_horvath)
horvath_coef <- horvath_coef_raw[, 1:2]

use_data(horvath_coef)
