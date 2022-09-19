# # Test methylation processing functions for pe hg38 wgbs data
# ##load_all(); test_file("tests/testthat/test-camdac_wgbs_methylation.R")
# 
# test_that("pipeline_files_load", {
# 
#   run_config <- create_camdac_config(
#     camdac_refs="./data/pipeline_files",
#     outdir="./runs/",
#     build="hg38",
#     bsseq="wgbs",
#     bsseq_lib="pe",
#     n_cores=5
#   )
#   
#   refs = get_reference_files(run_config, "loci_files")
#   testthat::expect_true(fs::file_exists(refs[[1]]))
#   
#   refs = get_reference_files(run_config, "segments_files")
#   testthat::expect_true(fs::file_exists(refs[[1]]))
#   
#   refs = get_reference_files(run_config, "annotations")
#   testthat::expect_true(fs::file_exists(refs[[1]]))
# }
# )
