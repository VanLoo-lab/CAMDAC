test_that("allele counting regions can be read from BED file", {
  # Est. 5 minutes
  bam <- system.file("testdata", "tumor.bam", package = "CAMDAC")
  seg_regions <- system.file("testdata", "test_segments.bed", package = "CAMDAC")

  # Create test config for segments only
  config_t <- config
  config_t$outdir <- "./result_test_segments"
  config_t$regions <- seg_regions

  # Run allele counting
  tumor <- CamSample(id = "T", sex = "XY", bam = bam)
  cmain_count_alleles(tumor, config)
  ac_file <- get_fpath(tumor, config, "counts")

  testthat::expect_true(
    fs::file_exists(
      ac_file
    )
  )
})
