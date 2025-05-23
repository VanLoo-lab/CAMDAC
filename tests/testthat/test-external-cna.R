# FIX:
test_that("CAMDAC deconvolves using pre-computed CNA data", {
  # Load external CNA object
  cna_file <- system.file("testdata", "test.cna.txt", package = "CAMDAC")

  # Delete expected CNA file if it exists
  exp <- get_fpath(tumor, config, "cna")
  if (fs::file_exists(exp)) {
    fs::file_delete(exp)
  }

  # Objects for tumor, normal and config defined in `tests/testthat/setup.R`
  # Add expected CNA file
  attach_output(tumor, config, "cna", cna_file)

  # Test that pipeline runs without calling CNAs
  stdout <- testthat::capture_output(
    pipeline(tumor, germline = NULL, infiltrates = normal, origin = normal, config)
  )

  testthat::expect_true(
    stringr::str_detect(stdout, "CNA Found.")
  )

  testthat::expect_true(
    stringr::str_detect(stdout, "pipeline complete")
  )
})
