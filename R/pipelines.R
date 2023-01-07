#' CAMDAC tumor-normal pipeline
#'
#' Run CAMDAC analysis on a bulk tumor and patient-matched tissue-matched tumor-adjacent normal sample.
#'
#' @param tumor. Tumor sample data built with `create_camdac_sample()`.
#' @param normal. Normal sample data built with `create_camdac_sample()`.
#' @param config. Configuration built with `create_camdac_config()`.
#' @export
pipeline_tumor_normal <- function(tumor, normal, config) {
  # Log
  loginfo("CAMDAC:::pipeline_tumor_normal start for %s", tumor$patient_id)

  # Count SNP and CpG alleles
  cmain_count_alleles(tumor, config)
  cmain_count_alleles(normal, config)

  # TODO: Refactor into a CMAIN function that accepts battenberg or ascat as input
  # If no allele-specific segments are given, call them
  # Or refit if cna object contains purity and ploidy values alone
  if (is.null(tumor$cna$ascna)) {
    # Call copy number abberations using battenberg
    cmain_make_snp_profiles(tumor, normal, config)
    cmain_run_battenberg(tumor, normal, config)
    tumor <- load_cna_data(tumor, config, "battenberg")
  } else {
    loginfo("Using external CNA data provided for %s", tumor$patient_id)
  }

  # Deconvolve methylation rates
  cmain_make_methylation_profile(tumor, config)
  cmain_make_methylation_profile(normal, config)
  cmain_deconvolve_methylation(tumor, normal, config)

  # Call differential methylation
  cmain_call_dmps(tumor, normal, config)
  cmain_call_dmrs(tumor, config)

  # Log
  loginfo("CAMDAC:::pipeline_tumor_normal complete for %s", tumor$patient_id)
}
