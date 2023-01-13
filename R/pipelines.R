#' CAMDAC tumor-normal pipeline
#'
#' Run CAMDAC analysis on a bulk tumor and patient-matched tissue-matched tumor-adjacent normal sample.
#'
#' @param tumor. Tumor sample data built with `CamSample()`.
#' @param normal. Normal sample data built with `CamSample()`.
#' @param config. Configuration built with `CamConfig()`.
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

#' CAMDAC tumor-panel pipeline
#'
#' Run CAMDAC analysis on a bulk tumor using a normal panel for deconvolution.
#' Limitations: Expects CNA data to be provided and methylation object built for normal.
#' @param tumor. Tumor sample data built with `create_camdac_sample()`. Requires attached CNA solution.
#' @param infiltrates. Normal sample data built with `create_camdac_sample()`. Must have methylation
#' object built or attached using `attach_methylation_panel()`.
#' @param cell_of_origin. Normal sample data built with `create_camdac_sample()`. Must have
#' methylation object attached using `attach_methylation_panel()`.
#' @param config. Configuration object created with `CamConfig()`.
#' @export
pipeline_tumor_panel <- function(tumor, infiltrates = normal, cell_of_origin = normal, config) {
  # Log
  loginfo("CAMDAC:::pipeline_tumor_panel start for %s", tumor$patient_id)

  # Perform validation checks
  validate_pipeline_tumor_panel(tumor, infiltrates, config)
  validate_pipeline_tumor_panel(tumor, cell_of_origin, config)

  # Run pipeline
  # Process tumor counts
  cmain_count_alleles(tumor, config)
  cmain_make_methylation_profile(tumor, config)

  # Deconvolve methylation rates
  cmain_deconvolve_methylation(tumor, infiltrates, config)

  # Call differential methylation
  cmain_call_dmps(tumor, cell_of_origin, config)
  cmain_call_dmrs(tumor, config)

  # Log
  loginfo("CAMDAC:::pipeline_tumor_panel complete for %s", tumor$patient_id)
}

validate_pipeline_tumor_panel <- function(tumor, normal, config) {
  # Check if tumor has CNA computed
  # TODO: Add tumor CNA attachment function as with methylation
  if (is.null(tumor$cna$ascna)) {
    stop("Tumor CNA data not found. Please add cna object to tumor `attach_cna_data()`.")
  }

  # Check if normal has methylation object
  normal_meth_exists <- file.exists(get_fpath(normal, config, "methylation"))
  if (!normal_meth_exists) {
    stop("Normal methylation data not found. Please run `attach_methylation_panel()` first.")
  }
}
