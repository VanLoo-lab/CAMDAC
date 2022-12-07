test_that("pipeline file download works", {
  pf <- download_pipeline_files("test", directory="pfdw_test")
  expected <- fs::path(pf, "wgbs")
  expect_true(fs::dir_exists(expected))
  fs::dir_delete(pf)
})
