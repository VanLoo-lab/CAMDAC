test_that("ASM pipeline runs", {
  # Run main

  main_config <- CamConfig(
    outdir = "./result_test",
    bsseq = "wgbs",
    build = "hg38",
    lib = "pe",
    regions = regions, # Speed up tests
    n_cores = 10,
    min_cov = 1, # Required to capture sufficient SNPs from test
    cna_caller = "ascat" # Battenberg always recommended. ASCAT used here to speed up test
  )

  # Run main (skips if outputs exist)
  pipeline(tumor, germline = normal, infiltrates = NULL, origin = NULL, config)

  asm_pipeline(
    tumor = tumor,
    germline = normal,
    infiltrates = normal,
    origin = normal,
    config = config
  )

  # Confirm that final output file is created
  asm_out <- get_fpath(tumor, asm_config, "asm_dmp")
  expect_true(file.exists(asm_out))
})
