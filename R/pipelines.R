#' CAMDAC tumor-normal pipeline
#'
#' Run CAMDAC analysis on a bulk tumor and patient-matched tissue-matched tumor-adjacent normal sample.
#'
#' @param tumor. Tumor sample data built with `CamSample()`.
#' @param normal. Normal sample data built with `CamSample()`.
#' @param config. Configuration built with `CamConfig()`.
#' @export
pipeline <- function(tumor, germline = NULL, infiltrates = NULL, origin = NULL, config) {
  # Log
  loginfo("CAMDAC:::pipeline start for %s", tumor$patient_id)

  # Preprocess CpG, SNP and methylation data for all samples
  preprocess(
    list(tumor, germline, infiltrates, origin),
    config
  )

  # Combine tumor-germline SNPs and call CNAs
  cmain_bind_snps(tumor, germline, config)
  cmain_call_cna(tumor, germline, config)

  # Run deconvolution
  cmain_deconvolve_methylation(tumor, infiltrates, config)

  # Call differential methylation
  cmain_call_dmps(tumor, origin, config)
  cmain_call_dmrs(tumor, config)

  # Log
  loginfo("CAMDAC:::pipeline complete for %s", tumor$patient_id)
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

#' Preprocess a list of CamSample objects for analysis
#' @param sample_list. List of CamSample objects.
#' @param config. CamConfig object.
#' @export
preprocess <- function(sample_list, config) {
  for (s in sample_list) {
    # Go to next part of loop if its null
    if (is.null(s)) {
      next
    }

    # Count SNP and CpG alleles if a BAM file is provided
    cmain_count_alleles(s, config)

    # Prepare SNP data for CNA calling if allele counts are present
    cmain_make_snps(s, config)

    # Format methylation rates for deconvolution
    cmain_make_methylation_profile(s, config)
  }
}
