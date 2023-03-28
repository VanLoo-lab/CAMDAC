test_that("allele counts combine to form panels given sample proportions", {
  # Test panel using two samples that should return 1 eligible CpG
  # Confirm that the panel methylation is a mixture of the two samples:
  #   sample 1 has a methylation of 1 and will be present at 20%
  #   sample 2 has a methylation of 0.5 and will be present at 80%
  #   expected pnael methylation at this site should be 0.6
  #   expected panel coverage for this site should be 20% from sample 1 and 80% from sample 2
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

# Test panel can be built from a matrix of beta values
data <- data.table::fread(
  system.file("testdata", "test_panel_from_beta.csv", package = "CAMDAC")
)
mat = data[, 4:ncol(data)]

panel_beta <- panel_meth_from_beta(
  mat = mat,
  chrom = data$chrom,
  start = data$start,
  end = data$end, 
  cov = 100, # Single value for coverage given to all CpGs.
  props = c(0.1, 0.8, 0.1),
  min_samples = 1,
  max_sd = 1
)

test_cg3 <- panel_beta[chrom == "13" & start == 18231437 & end == 18231438, ]
expect_equal(round(test_cg3$m, 2), 0.24) # Expect linear combination of three betas
expect_equal(test_cg3$M + test_cg3$UM, 100) # Expect cg cov to meet input value
test_cg4 <- panel_beta[chrom == "13" & start == 18173666 & end == 18173667, ]
expect_equal(round(test_cg4$m, 2), 0.15)

# Test panel can be created from a matrix of beta values and a vector 
data <- data.table::fread(
  system.file("testdata", "test_panel_from_beta.csv", package = "CAMDAC")
)[1:5,]
mat = data[, 4:ncol(data)]

panel_cov <- panel_meth_from_beta(
  mat = mat,
  chrom = data$chrom,
  start = data$start,
  end = data$end, 
  cov = c(10, 1, 10, 1, 10), # Single value for coverage given to all CpGs.
  props = c(0.1, 0.8, 0.1),
  min_samples = 1,
  max_sd = 1
)

panel_cov2 <- panel_meth_from_beta(
  mat = mat,
  chrom = data$chrom,
  start = data$start,
  end = data$end, 
  cov = matrix(30, nrow=5, ncol=3), # Single value for coverage given to all CpGs.
  props = c(0.1, 0.8, 0.1),
  min_samples = 1,
  max_sd = 1
)

expect_equal(panel_cov$cov, c(10, 1, 10, 1, 10))
expect_equal(panel_cov2$cov, c(30, 30, 30, 30, 30))

})
