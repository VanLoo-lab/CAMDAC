test_that("single-end allele counter runs on hg19 samples samples", {
    # CNA caller test config
    config_hg19 <- CamConfig(
        outdir = "./result_hg19",
        bsseq = "wgbs",
        build = "hg19",
        lib = "se",
        n_cores = 10,
        min_cov = 1 # Required to capture sufficient SNPs from test
    )

    # Test that beagle found for hg19
    testthat::expect_true(file.exists(config_hg19$beaglejar))
    # Test that caller is not battenberg (not implemented yet for hg19)
    testthat::expect_true(config_hg19$cna_caller != "battenberg")

    # Run allele count on a small segment and confirm that BAM file with non `chr` mapping is processed
    # bam <- system.file("testdata", "normal_hg19.bam", package = "CAMDAC")
    # seg_regions <- system.file("testdata", "test_segments.bed", package = "CAMDAC")
    # normal <- CamSample(id = "N", sex = "XY", bam = bam)
    # ac_file <- get_fpath(normal, config_hg19, "counts")
    # config_hg19$regions <- seg_regions
    # cmain_count_alleles(normal, config_hg19)
    # testthat::expect_true(fs::file_exists(ac_file))
    # fs::file_delete(ac_file) # Remove file to use new segments without overwriting.

    # Test full pipeline on hg19 
    ## Overwrite regions for full test
    regions <- system.file("testdata", "test_wgbs_segments.bed", package = "CAMDAC")
    config_hg19$regions = regions

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
