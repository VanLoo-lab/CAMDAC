test_that("Pipeline tumour-normal completes", {
    testthat::skip("Skipped. Long-running.")

    bam <- system.file("testdata", "tumor.bam", package = "CAMDAC")
    bam2 <- system.file("testdata", "normal.bam", package = "CAMDAC")
    regions <- system.file("testdata", "test_wgbs_small.bed", package = "CAMDAC")

    config <- CamConfig(
        outdir = "./result_test",
        bsseq = "wgbs",
        build = "hg38",
        lib = "pe",
        regions = regions,
        n_cores = 10
    )

    tumor <- CamSample(id = "T", sex = "XY", bam = bam)
    normal <- CamSample(id = "N", sex = "XY", bam = bam2)

    stdout <- testthat::capture_output(
        pipeline(tumor, germline = normal, infiltrates = normal, origin = normal, config)
    )

    expr <- "pipeline_tumor_normal complete"

    testthat::expect_true(
        stringr::str_detect(stdout, expr)
    )
})
