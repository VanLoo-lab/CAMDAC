test_that("ascat and battenberg runs on wgbs samples", {
    # Preprocess CpG, SNP and methylation data for all samples
    preprocess(
        list(tumor, normal),
        config
    )

    # Combine tumor-germline SNPs and call CNAs
    cmain_bind_snps(tumor, normal, config)

    # Test CNA and expect file exists after ASCAT
    cna_file <- get_fpath(tumor, config, "cna")
    if (fs::file_exists(cna_file)) {
        fs::file_delete(cna_file)
    }
    config$cna_caller <- "ascat"
    cmain_call_cna(tumor, normal, config)

    # Reset
    testthat::expect_true(fs::file_exists(cna_file))
    fs::file_delete(cna_file)

    # Run battenberg
    config$cna_caller <- "battenberg"
    # Battenberg warnings are function of the test data
    suppressWarnings(cmain_call_cna(tumor, normal, config))
    testthat::expect_true(fs::file_exists(cna_file))
})
