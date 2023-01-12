test_that("CAMDAC deconvolves using pre-computed CNA data", {
  # Load external CNA object
  cna <- readRDS(system.file("testdata", "test_cna.rds", package = "CAMDAC"))
  segments_bed <- system.file("testdata", "test_wgbs_small_segments.bed", package = "CAMDAC")

  # Warning: Currently requires 30 minutes to run and CAMDAC_PIPELINE_FILES environment var
  #    to be set to the location of the CAMDAC pipeline files.
  config <- CamConfig(
    outdir = "./result_test", # Path for outputs
    bsseq = "wgbs", # WGBS
    build = "hg38", # Reference
    bsseq_lib = "pe", # Paired end
    n_cores = 10 # Cores for parallel processing
  )

  # Get test BAM
  tumor_bam <- system.file("extdata", "NA18939_test_v3.bam", package = "CAMDAC")
  normal_bam <- system.file("extdata", "NA20502_test_v3.bam", package = "CAMDAC")

  # Set CAMDAC tumor sample
  tumor <- CamSample(
    patient_id = "P1",
    patient_sex = "XX",
    sample_id = "T",
    sample_type = "tumour",
    bam_file = tumor_bam,
    segments_bed = segments_bed,
    cna = cna
  )

  # Set CAMDAC normal sample
  normal <- CamSample(
    patient_id = "P1",
    patient_sex = "XX",
    sample_id = "N",
    sample_type = "normal",
    bam_file = normal_bam,
    segments_bed = segments_bed
  )

  # Run pipeline
  stdout <- testthat::capture_output(
    pipeline_tumor_normal(tumor, normal, config)
  )

  expr_cna <- "Using external CNA data"
  expr <- "pipeline_tumor_normal complete"

  testthat::expect_true(
    stringr::str_detect(stdout, expr_cna)
  )

  testthat::expect_true(
    stringr::str_detect(stdout, expr)
  )
})
