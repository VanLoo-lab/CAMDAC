test_that("tumor panel pipeline runs with panel of normals", {
  # Create test config
  config_p <- config
  config_p$outdir <- "result_panel"

  # Get internal datasets
  cna_file <- system.file("testdata", "test.cna.txt", package = "CAMDAC")
  meth_file <- system.file("testdata", "test_panel.m.csv.gz", package = "CAMDAC")

  # Setup panel
  panel <- CamSample(id = "PANEL", sex = "XY", bam = NULL)
  attach_output(panel, config_p, "meth", meth_file)
  attach_output(tumor, config_p, "cna", cna_file)

  # Test pipeline runs
  stdout <- testthat::capture_output(
    pipeline(tumor, germline = normal, infiltrates = panel, origin = panel, config_p)
  )

  expr <- "pipeline complete"
  testthat::expect_true(
    stringr::str_detect(stdout, expr)
  )
})
