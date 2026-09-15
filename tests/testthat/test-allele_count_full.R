test_that("Tumor allele counting completes for test sample", {
  bam <- system.file("testdata", "tumor.bam", package = "CAMDAC")

  # Create test config for segments only
  config_t <- config
  config_t$outdir <- "./result_test_ac_full"
  withr::defer(fs::dir_delete(config_t$outdir))

  # Run allele counting
  tumor <- CamSample(id = "T", sex = "XY", bam = bam)
  cmain_count_alleles(tumor, config_t)
  ac_file <- get_fpath(tumor, config_t, "counts")

  testthat::expect_true(
    fs::file_exists(
      ac_file
    )
  )
})
