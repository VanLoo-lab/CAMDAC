test_that("tumor panel pipeline runs with panel of normals", {
  # Get internal datasets
  cna <- readRDS(system.file("testdata", "test_cna.rds", package = "CAMDAC"))
  segments_bed <- system.file("testdata", "test_wgbs_small_segments.bed", package = "CAMDAC")
  tumor_bam <- system.file("extdata", "NA18939_test_v3.bam", package = "CAMDAC")

  config <- create_camdac_config(
    outdir = "./panel_test", # Path for outputs
    bsseq = "wgbs", # WGBS
    build = "hg38", # Reference
    bsseq_lib = "pe", # Paired end
    n_cores = 10 # Cores for parallel processing
  )

  # Load Panel
  panel <- fread(system.file("testdata", "test_panel.csv.gz", package = "CAMDAC"))

  # Load tumor from BAM with pre-computed CNAs
  tumor <- create_camdac_sample(
    patient_id = "P1",
    patient_sex = "XX",
    sample_id = "T",
    sample_type = "tumour",
    bam_file = tumor_bam,
    segments_bed = segments_bed,
    cna = cna
  )

  # Load normal from panels
  normal <- create_camdac_sample(
    patient_id = "P1",
    patient_sex = "XX",
    sample_id = "PANEL",
    sample_type = "normal",
    bam_file = NULL # No BAM given for panel object
  )
  normal <- attach_methylation_panel(normal, config, panel)

  stdout <- testthat::capture_output(
    pipeline_tumor_panel(tumor, infiltrates = normal, cell_of_origin = normal, config)
  )

  expr <- "pipeline_tumor_panel complete"

  testthat::expect_true(
    stringr::str_detect(stdout, expr)
  )
})
