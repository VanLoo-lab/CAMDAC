
#' Make CAMDAC methylation panel from allele counts
#' @param ac_files Allele count files from CAMDAC
#' @param min_coverage Minimum coverage for a sample's site to be included in panel
#' @param min_samples Minimum number of samples with coverage for a site to be included in panel
#' @param max_sd Maximum standard deviation of methylation for a site to be included in panel
#' @param drop_snps Boolean. If TRUE, drop per-sample CG-SNPs (BAF < 0.1 or BAF > 0.9) from panel
#' @export
panel_meth_from_counts <- function(ac_files, min_coverage = 3, min_samples = 1,
                                   max_sd = 0.1, drop_snps = FALSE) {
  # Load AC files

  # Find unique CpG positions

  # Order AC files by unique CG overlap

  # Filter CpGs by minimum sample coverage threshold

  # Mask sample CpGs by coverage

  # Mask sample CpGs by SNP loci if required

  # Filter CpGs by maximum SD threshold

  # Combine counts to create methylation panel

  # Return panel object

  # Test pipeline:
  # load_all();test_active_file("tests/testthat/test-panel.R")
  x <- data.table(total_counts_m = 100)
  x[, `:=`("chrom" = 1, "start" = 2, "end" = 3, "M_n" = 10, "UM_n" = 10, "m_n" = 10, "cov_n" = 5)]
  return(x)
}
