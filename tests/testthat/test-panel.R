test_that("allele counts combine to form panels", {
  
  # load test allele counts
  ac_file <- system.file("extdata", "NA18939_test_ac.SNPs.CpGs.all.sorted.csv.gz", package="CAMDAC")
  panel <- panel_meth_from_counts(
    ac_files=c(ac_file, ac_file, ac_file),
    min_coverage=3,
    min_samples=1,
    max_sd=0.1,
    drop_snps=T
    )
  
  # Test panel is a data table
  expect_is(panel, "data.table")
  # Test panel has expected rows
  expect_true(nrow(panel) >= 1)
  # Test panel has the expected columns
  expected_fields <- c("chrom", "start", "end", "M_n", "UM_n", "m_n", "cov_n")
  expect_true(
    all(expected_fields %in% names(panel))
  )

  # Test panel has more reads than an individual file
  ac <- data.table::fread(ac_file)
  counts_bool = (ac$total_counts_m < panel$total_counts_m) | (is.na(ac$total_counts_m))
  expect_true(all(counts_bool))

})
