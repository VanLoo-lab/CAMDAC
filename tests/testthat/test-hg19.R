test_that("CAMDAC runs with battenberg in hg19 mode", {
    # CNA caller test config
    config_hg19 <- CamConfig(
        outdir = "./result_hg19_pe",
        bsseq = "wgbs",
        build = "hg19",
        lib = "pe",
        cna_caller = "battenberg",
        n_cores = 10,
        min_cov = 1 # Required to capture sufficient SNPs from test
    )
    config_hg19$regions = system.file("testdata", "test_wgbs_segments.bed", package = "CAMDAC")

    stdout <- testthat::capture_output(
        pipeline(tumor, germline = normal, infiltrates = normal, origin = normal, config_hg19)
    )

    expr <- "pipeline complete"

    testthat::expect_true(
        stringr::str_detect(stdout, expr)
    )

})
