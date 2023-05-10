test_that("ascat and battenberg runs on wgbs samples", {
    # CNA caller test config
    config_c <- CamConfig(
        outdir = "./result_cna",
        bsseq = "wgbs",
        build = "hg38",
        lib = "pe",
        regions = regions,
        n_cores = 30,
        min_cov = 1 # Required to capture sufficient SNPs from test
    )
    withr::defer(fs::dir_delete(config_c$outdir))

    # Attach tsnps file to tumor
    tsnps_file <- system.file("testdata", "test_tsnps.csv.gz", package = "CAMDAC")
    attach_output(tumor, config_c, "tsnps", tsnps_file)

    # Reset
    cna_file <- get_fpath(tumor, config_c, "cna")
    if (fs::file_exists(cna_file)) {
        fs::file_delete(cna_file)
    }

    # Test CNA and expect file exists after ASCAT
    config_c$cna_caller <- "ascat"
    cmain_call_cna(tumor, normal, config_c)

    # Test ASCAT
    tool <- fread(cna_file)$pipeline[[1]]
    testthat::expect_equal(tool, "ascat")
    fs::file_delete(cna_file)

    # Run battenberg
    config_c$cna_caller <- "battenberg"
    cmain_count_alleles(tumor, config_c)
    cmain_count_alleles(normal, config_c)

    # Battenberg warnings are function of the test data
    suppressWarnings(cmain_call_cna(tumor, normal, config_c))
    tool <- fread(cna_file)$pipeline[[1]]
    testthat::expect_equal(tool, "battenberg")
    fs::file_delete(cna_file)
})
