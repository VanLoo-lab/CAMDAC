test_that("Pipeline tumour-normal completes", {
    testthat::skip("Skipped. Long-running.")

    # Warning: Currently requires 30 minutes to run and CAMDAC_PIPELINE_FILES environment var
    #    to be set to the location of the CAMDAC pipeline files.
    config <- CamConfig(
        outdir = "./result_test", # Path for outputs
        bsseq = "wgbs", # WGBS
        build = "hg38", # Reference
        bsseq_lib = "pe", # Paired end
        n_cores = 10 # Cores for parallel processing
    )

    # Get test BAM
    tumor_bam <- system.file("extdata", "NA18939_test_v3.bam", package = "CAMDAC")
    normal_bam <- system.file("extdata", "NA20502_test_v3.bam", package = "CAMDAC")

    # Set CAMDAC tumor sample
    tumor <- CamSample(
        patient_id = "P1",
        patient_sex = "XX",
        sample_id = "T",
        sample_type = "tumour",
        bam_file = tumor_bam
    )

    # Set CAMDAC normal sample
    normal <- CamSample(
        patient_id = "P1",
        patient_sex = "XX",
        sample_id = "N",
        sample_type = "normal",
        bam_file = normal_bam
    )

    stdout <- testthat::capture_output(
        pipeline_tumor_normal(tumor, normal, config)
    )

    expr <- "pipeline_tumor_normal complete"

    testthat::expect_true(
        stringr::str_detect(stdout, expr)
    )
})
