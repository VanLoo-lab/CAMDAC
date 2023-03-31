test_that("ASM pipeline runs", {
  # Setup config
  asm_config <- CamConfig(
    outdir = "./result_asm_full", bsseq = "wgbs", lib = "pe",
    build = "hg38", n_cores = 3, min_cov = 1, cna_caller = "ascat"
  )

  # Add ASM CNA caller
  attach_output(tumor, asm_config, "asm_cna", system.file("testdata", "tumor.cna.txt", package = "CAMDAC"))
  attach_output(tumor, asm_config, "asm_snps", system.file("testdata", "test_het_snps.tsv", package = "CAMDAC"))
  attach_output(normal, asm_config, "asm_snps", system.file("testdata", "test_het_snps.tsv", package = "CAMDAC"))

  asm_pipeline(
    tumor = tumor,
    germline = normal,
    infiltrates = normal,
    origin = normal,
    config = asm_config
  )

  # Confirm that final output file is created
  asm_out <- get_fpath(tumor, asm_config, "asm_dmp")
  expect_true(file.exists(asm_out))
})
