test_that("ASM allele counter runs", {
  # Load hets
  # Attach SNPs to tumor and normal objects
  hets_file <- system.file("testdata", "test_het_snps.tsv", package="CAMDAC")

  asm_config <- CamConfig(
    outdir = "./result_asm",
    bsseq = "wgbs",
    build = "hg38",
    lib = "pe",
    n_cores = 10,
    min_cov = 1 # Required to capture sufficient SNPs from test
  )

  attach_output(tumor, asm_config, "asm_snps", hets_file)
  attach_output(normal, asm_config, "asm_snps", hets_file)

  # Run ASM allele counter
  cmain_asm_allele_counts(tumor, asm_config)
  cmain_asm_allele_counts(normal, asm_config)

  # Test for expected output files
  exp_files <- c(
    get_fpath(tumor, asm_config, "asm_counts"),
    get_fpath(tumor, asm_config, "asm_hap_stats"),
    get_fpath(tumor, asm_config, "asm_phase_map"),
    get_fpath(normal, asm_config, "asm_counts"),
    get_fpath(normal, asm_config, "asm_hap_stats"),
    get_fpath(normal, asm_config, "asm_phase_map")
  )
  testthat::expect_true(all(file.exists(exp_files)))
})
