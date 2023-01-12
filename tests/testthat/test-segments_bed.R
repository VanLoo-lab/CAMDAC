test_that("allele counting regions can be read from BED file", {
  # Load test sample and config
  config <- create_camdac_config(
    # camdac_refs = "./camdac_refs", # Optional: Pass the location of pipeline files
    outdir = "./wgbs_test", # Path for outputs
    bsseq = "wgbs", # WGBS
    build = "hg38", # Reference
    bsseq_lib = "pe", # Paired end
    n_cores = 10 # Cores for parallel processing
  )

  tumor <- CamSample(
    patient_id = "P1",
    patient_sex = "XX",
    sample_id = "T",
    sample_type = "tumour",
    # Load segments from BED file
    bam_file = system.file("extdata", "NA18939_test_v3.bam", package = "CAMDAC"),
    segments_bed = system.file("testdata", "test_segments.bed", package = "CAMDAC")
  )

  cmain_count_alleles(tumor, config)
  ac_file <- build_output_name(tumor, config, "allele_counts")

  testthat::expect_true(
    fs::file_exists(
      ac_file
    )
  )
})
