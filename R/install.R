# List of packages required
list.of.packages <- c("optparse", "stringr", "GenomicRanges", "Rsamtools", "scales", 
                      "ggplot2", "gridExtra", "gtable", "data.table","readr", "dplyr", 
                      "parallel", "doParallel", "devtools", "fst")

# Check for bioconductors
if(!"BiocManager"%in% installed.packages()[,"Package"]) install.packages("BiocManager")
library(BiocManager)
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) BiocManager::install(new.packages)
 
# Check for CRAN
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

# Load ASCAT lib
library(devtools)
devtools::install_github("VanLoo-lab/ascat/ASCAT@v2.5.3")
#devtools::install_github("fstpackage/fst")
