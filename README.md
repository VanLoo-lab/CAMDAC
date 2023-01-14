
<!-- README.md is generated from README.Rmd. Please edit that file -->

# CAMDAC

Copy-number Aware Methylation Deconvolution Analysis of Cancer (CAMDAC)
is an R library for deconvolving pure tumor methylation from bulk tumor
sequencing data ([Larose Cadieux et al., 2022,
bioXriv](https://www.biorxiv.org/content/10.1101/2020.11.03.366252v2)).

This branch performs CAMDAC analysis for whole genome bisulfite
sequencing (WGBS) data. For Reduced Representation Bisulfite Sequencing
(RRBS) data analysis, visit the [VanLoo-lab/CAMDAC main
branch](https://github.com/VanLoo-lab/CAMDAC/tree/main).

<!-- badges: start -->

<!-- badges: end -->

## Documentation

View the full documentation at <https://vanloo-lab.github.io/CAMDAC/>.

## Quickstart

CAMDAC can be installed from an R console:

``` r
remotes::install_github("VanLoo-lab/CAMDAC@wgbs")
CAMDAC::download_pipeline_files("wgbs")
```

To run the tumor-normal deconvolution pipeline with test data:

``` r
library(CAMDAC)

tumor_bam = system.file("testdata", "tumor.bam", package = "CAMDAC")
normal_bam = system.file("testdata", "normal.bam", package = "CAMDAC")

tumor = CamSample(id="T", sex="XY", bam=tumor_bam)
normal = CamSample(id="N", sex="XY", bam=normal_bam)
config = CamConfig(outdir="./results", bsseq="wgbs", lib="pe", build="hg38")
CAMDAC::pipeline(tumor, germline=normal, infiltrates=normal, origin=normal, config)
```

# Development

``` r
library(devtools)
devtools::install_dev_deps("VanLoo-lab/CAMDAC@wgbs")
# make 
```
