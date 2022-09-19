
test_that("Methylation profile is created", {
  tumor <- CAMDAC::create_camdac_sample(
    patient_id="TEST",
    patient_sex="XX",
    sample_id="T",
    sample_type="tumour",
    bam_file=system.file("extdata", "NA18939_test_v3.bam", package = "CAMDAC")
  )
  
  test_result = system.file("extdata", "test_result", package = "CAMDAC")
  config <- CAMDAC::create_camdac_config(
    camdac_refs=".",
    outdir=test_result,
    build="hg38",
    bsseq="wgbs",
    bsseq_lib="pe"
  )
  
  expected_output = CAMDAC::build_output_name(tumor, config, "methylation")
  
  CAMDAC::cmain_make_methylation_profile(tumor, config)
  testthat::expect_true(
    fs::file_exists(
      expected_output
    )
  )
})
