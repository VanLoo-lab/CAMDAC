# Global text fixtures that may be used by tests.
bam <- system.file("testdata", "tumour_beds_min.sorted.bam", package = "CAMDAC")
bam2 <- system.file("testdata", "normal_beds_min.sorted.bam", package = "CAMDAC")
regions <- system.file("testdata", "test_beds_segments.bed", package = "CAMDAC")

config <- CamConfig(
    outdir = "./result_test",
    bsseq = "wgbs",
    build = "hg38",
    lib = "pe",
    regions = regions, # Speed up tests
    n_cores = 10,
    min_cov = 1, # Required to capture sufficient SNPs from test
    min_normal_cov = 1,
    min_mapq = 1
)

tumor <- CamSample(id = "T", sex = "XY", bam = bam)
normal <- CamSample(id = "N", sex = "XY", bam = bam2)

# Cleanup, as presented in https://testthat.r-lib.org/articles/test-fixtures.html
withr::defer(teardown_env())
