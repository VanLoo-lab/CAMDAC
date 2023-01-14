test_that("attach output function writes files", {
  outdir <- tempdir()

  config <- CamConfig(
    outdir = outdir,
    bsseq = "wgbs",
    build = "hg38",
    lib = "pe",
    n_cores = 10
  )

  tumor <- CamSample(id = "T", sex = "XY", bam = bam)

  # Get expected file and test that it doesn't exist
  exp <- get_fpath(tumor, config, "counts")
  if (file.exists(exp)) {
    fs::file_delete(exp)
  }

  # Run attach_output
  counts_file <- system.file("testdata", "test.SNPs.CpGs.all.sorted.csv.gz", package = "CAMDAC")
  attach_output(tumor, config, "counts", counts_file)

  # Test that file exists
  testthat::expect_true(file.exists(exp))
})
