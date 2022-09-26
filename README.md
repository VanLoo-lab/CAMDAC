
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CAMDAC

Copy-number Aware Methylation Deconvolution Analysis of Cancer (CAMDAC)
is a tool that deconvolves pure tumor methylation from bulk tumor
sequencing data ([Larose Cadieux et al., 2022,
bioXriv](https://www.biorxiv.org/content/10.1101/2020.11.03.366252v2)).

This branch describes CAMDAC for whole genome bisulfite sequencing
(WGBS) data. To run CAMDAC on Reduced Representation Bisulfite
Sequencing (RRBS) data, visit the [VanLoo-lab/CAMDAC main
branch](https://github.com/VanLoo-lab/CAMDAC/tree/main).

<!-- badges: start -->

<!-- badges: end -->

## Installation

CAMDAC can be installed from an R console:

``` r
remotes::install_github("VanLoo-lab/CAMDAC@wgbs")
CAMDAC::download_pipeline_files("wgbs")
```

Additional dependencies:

  - java runtime environment [(jre)](https://openjdk.org/)
  - [beagle5](http://faculty.washington.edu/browning/beagle/beagle.html)
    (included)

## Documentation

View the full documentation at <https://vanloo-lab.github.io/CAMDAC/>.

## Quickstart

CAMDAC tumor-normal deconvolution pipeline with test data:

``` r
library(CAMDAC)
tumor = create_camdac_sample("P1_T1", bam_file = "tumor.bam")
normal = create_camdac_sample("P1_N1", bam_file = "normal.bam")
config = create_camdac_config(outdir="./results", bsseq="wgbs", bsseq_lib="pe", build="hg38")
pipeline_tumor_normal(tumor, normal, config)
```
