# Copy-number Aware Methylation Deconvolution and Analysis of Cancer (CAMDAC)

Plesae refer to the [CAMDAC manual](https://htmlpreview.github.io/?https://github.com/VanLoo-lab/CAMDAC/blob/main/CAMDAC_manual/CAMDAC_manual.html) for a detailed description of the CAMDAC principles, installation and steps for running the code.

To cite CAMDAC, please refer to our pre-print: [Larose Cadieux et al., 2020. Copy number-aware deconvolution of tumor-normal DNA methylation profiles. bioRxiv.](https://doi.org/10.1101/2020.11.03.366252).

## Installation

The CAMDAC R library can be install from github repository:

```r
# Install the remotes package 
install.packages("remotes")

# Install CAMDAC from GitHub
remotes::install_github("VanLoo-lab/CAMDAC")
```

Files required to run the CAMDAC pipeline [(listed here)](inst/extdata/pipeline_files_urls.txt) can be downloaded with a helper function:

```r
library(CAMDAC)
CAMDAC::download_pipeline_files(bsseq="rrbs", directory="pipeline_files/")
```

## Quickstart

To call CAMDAC with a matched tumor and adjacent normal sample:

```r
library(CAMDAC)
tumor_bam <- system.file("extdata", "test_tumor.bam", package = "CAMDAC")
normal_bam <- system.file("extdata", "test_normal.bam", package = "CAMDAC")

CAMDAC::pipeline_tumor_normal(
    patient_id="P1",
    tumor_id="T",
    normal_id="N",
    tumor_bam=tumor_bam,
    normal_bam=normal_bam,
    sex="XY",
    path="results/",
    pipeline_files="pipeline_files/",
    build="hg38",
    min_tumor = 1,
    min_normal = 1,
    mq = 0,
    n_cores = 1
)
```
