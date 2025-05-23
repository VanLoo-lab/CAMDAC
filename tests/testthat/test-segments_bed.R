# FIX
test_that("allele counting regions can be read from BED file", {
  bam <- system.file("testdata", "tumor.bam", package = "CAMDAC")
  seg_regions <- system.file("testdata", "test_segments.bed", package = "CAMDAC")

  # Create test config for segments only
  config_t <- config
  config_t$outdir <- "./result_test_segments"
  config_t$regions <- seg_regions
  config_t$overwrite <- TRUE
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

  # Check that allele counts regions do not extend beyond the segments given
  dt <- fread_chrom(ac_file)
  regs <- data.table::fread(seg_regions)
  names(regs) <- c("chrom", "start", "end")

  testthat::expect_equal(
    unique(dt$CHR),
    unique(regs$chrom)
  )

  all_within <- max(dt$start) <= max(regs$end) && min(dt$end) >= min(regs$start)
  testthat::expect_true(
    all_within
  )
})
