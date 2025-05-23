# Global text fixtures that may be used by tests.
bam <- system.file("testdata", "tumor.bam", package = "CAMDAC")
bam2 <- system.file("testdata", "normal.bam", package = "CAMDAC")
regions <- system.file("testdata", "test_wgbs_segments.bed", package = "CAMDAC")

config <- CamConfig(
    outdir = "./result_test",
    bsseq = "wgbs",
    build = "hg38",
    lib = "pe",
    regions = regions, # Speed up tests
    n_cores = 10,
    min_cov = 1 # Required to capture sufficient SNPs from test
)

tumor <- CamSample(id = "T", sex = "XY", bam = bam)
normal <- CamSample(id = "N", sex = "XY", bam = bam2)

# Cleanup, as presented in https://testthat.r-lib.org/articles/test-fixtures.html
withr::defer(teardown_env())
