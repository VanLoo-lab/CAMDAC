test_that("ascat and battenberg runs on wgbs samples", {
    # Est. 5 minutes
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

    # Preprocess
    cmain_count_alleles(tumor, config)
    cmain_count_alleles(normal, config)
    cmain_make_snps(tumor, config)
    cmain_make_snps(normal, config)
    cmain_bind_snps(tumor, normal, config)

    # Test CNA and expect file exists after ASCAT
    cna_file <- get_fpath(tumor, config, "cna")
    cmain_run_ascat(tumor, normal, config)

    # Reset
    testthat::expect_true(fs::file_exists(cna_file))
    fs::file_delete(cna_file)
    testthat::expect_false(fs::file_exists(cna_file))

    # Run battenberg
    cmain_run_battenberg(tumor, normal, config)
    testthat::expect_true(fs::file_exists(cna_file))
})
