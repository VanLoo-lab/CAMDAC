test_that("CAMDAC runs on hg19 samples samples", {
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

    # Test full pipeline on hg19 
    ## N.B. : Regions only returns 19 sites. Need to rebuild hg19 test BAM
    ## Overwrite regions for full test
    ## regions <- system.file("testdata", "test_wgbs_segments_hg19_v2.bed", package = "CAMDAC")
    ## config_hg19$regions = regions

    tumor <- CamSample(id = "T", sex = "XY",
        bam = system.file("testdata", "tumor_hg19.bam", package = "CAMDAC"))
    normal <- CamSample(id = "N", sex = "XY",
        bam = system.file("testdata", "normal_hg19.bam", package = "CAMDAC"))  
    stdout <- testthat::capture_output(
        pipeline(tumor, germline = normal, infiltrates = normal, origin = normal, config_hg19)
    )

    expr <- "pipeline complete"

    testthat::expect_true(
        stringr::str_detect(stdout, expr)
    )

})
