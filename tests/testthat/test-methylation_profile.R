
test_that("Methylation profile is created", {
  tumor <- create_camdac_sample(
    patient_id="TEST",
    patient_sex="XX",
    sample_id="T",
    sample_type="tumour",
    bam_file=system.file("extdata", "NA18939_test_v3.bam", package = "CAMDAC")
  )
  
  #test_result = system.file("extdata", "test_result", package = "CAMDAC")
  test_result = "test_wgbs_out"
  config <- create_camdac_config(
    camdac_refs=".",
    outdir=test_result,
    build="hg38",
    bsseq="wgbs",
    bsseq_lib="pe"
  )
  
  #expected_output = build_output_name(tumor, config, "methylation")
  #cmain_make_methylation_profile(tumor, config)
  testthat::expect_true(
    fs::dir_exists(
      test_result
    )
  )
})
