test_that("allele counts combine to form panels given sample proportions", {
  # Test panel using two samples that should return 1 eligible CpG
  # Confirm that the panel methylation is a mixture of the two samples:
  #   sample 1 has a methylation of 1 and will be present at 20%
  #   sample 2 has a methylation of 0.5 and will be present at 80%
  #   expected pnael methylation at this site should be 0.6
  #   expected panel coverage for this site should be 20% reads from sample 1 and 80$ counts from sample 2
  ac_sample1 <- system.file("testdata", "test.SNPs.CpGs.all.sorted.csv.gz", package = "CAMDAC")
  ac_sample2 <- system.file("testdata", "test_prop.SNPs.CpGs.all.sorted.csv.gz", package = "CAMDAC")

  panel <- panel_meth_from_counts(
    ac_files = c(ac_sample1, ac_sample2),
    ac_props = c(0.2, 0.8),
    min_coverage = 3,
    min_samples = 1,
    max_sd = 0.8,
    drop_snps = TRUE
  )

  # Test panel is a data table
  expect_is(panel, "data.table")
  # Test panel has expected rows
  expect_true(nrow(panel) >= 1)
  # Test panel has the expected columns
  expected_fields <- c("chrom", "start", "end", "M", "UM", "m", "cov")
  expect_true(
    all(expected_fields %in% names(panel))
  )

  # Test panel has the expected methylation
  test_cg <- panel[chrom == "13" & start == 18231437 & end == 18231438, ]
  expect_equal(test_cg$m, 0.6)

  # Test panel has the expected coverage
  expect_equal(test_cg$cov, 12)

  # Test panel default behaviour is to return the sum of counts
  panel_default <- panel_meth_from_counts(
    ac_files = c(ac_sample1, ac_sample2),
    min_coverage = 3,
    min_samples = 1,
    max_sd = 0.9,
    drop_snps = TRUE
  )
  test_cg2 <- panel_default[chrom == "13" & start == 18231437 & end == 18231438, ]
  expect_equal(test_cg2$M, 8)
  expect_equal(test_cg2$UM, 4)
  expect_equal(test_cg2$m, (8 / (8 + 4)))
})
