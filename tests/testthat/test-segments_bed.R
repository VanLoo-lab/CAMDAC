test_that("allele counting regions can be read from BED file", {
  # Est. 5 minutes
  bam <- system.file("testdata", "tumor.bam", package = "CAMDAC")
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
  cmain_count_alleles(tumor, config)
  ac_file <- get_fpath(tumor, config, "counts")

  testthat::expect_true(
    fs::file_exists(
      ac_file
    )
  )
})
