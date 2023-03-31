#' CAMDAC tumor-normal pipeline
#'
#' Run CAMDAC analysis on a bulk tumor and patient-matched tissue-matched tumor-adjacent normal sample.
#'
#' @param tumor Tumor `CamSample` object for deconvultion.
#' @param germline Patient-matched normal `CamSample` object. May be NULL if `tumor` has CNA calls already.
#' @param infiltrates Normal `CamSample` as a proxy for infiltrating normal methylation.
#' @param origin Normal `CamSample` representing cell of origin for tumor-normal differential methylation.
#' @param config Configuration built with `CamConfig()`.
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

#' Run allele-specific methylation analysis pipeline
#' @param tumor. CamSample object for tumor sample.
#' @param germline. CamSample object for germline sample. Used for CNA calling.
#' @param infiltrates. CamSample object for infiltrating normal sample. Used for deconvolution.
#' @param origin. CamSample object for cell of origin sample. Used for differential methylation.
#' @param config. CamConfig object.
asm_pipeline <- function(tumor, germline = NULL, infiltrates = NULL, origin = NULL, config) {
  # Log
  loginfo("CAMDAC:::asm_pipeline start for %s", tumor$patient_id)
  sample_list <- list(tumor, germline, infiltrates, origin)

  # Checks that ASM SNPs file is available, otherwise, creates from bulk allele counts on germline
  cmain_asm_make_snps(tumor, germline, config)

  # Check that ASM CNA file is available, otherwise, create from CAMDAC CNA calls
  cmain_asm_call_cna(tumor, germline, config)

  # Preprocess CpG, SNP and methylation data for all samples
  loginfo("Preprocessing ASM data")
  preprocess_asm(
    sample_list,
    config
  )

  # Assign ASM CNA to per-allele CG sites
  cmain_fit_meth_cna(tumor, config)

  # Run ASM deconvolution
  cmain_asm_deconvolve(tumor, infiltrates, config)

  # TODO: How are CG-SNPs handled at allele counting stage for ASM?
  # Run ASM differential methylation within-sample
  for (s in sample_list) {
    if (!is.null(s)) {
      cmain_asm_ss_dmps(s, config)
    }
  }

  # Run ASM differential methylation between samples
  cmain_asm_dmps(tumor, origin, config)

  # Log complete
  loginfo("CAMDAC:::asm_pipeline complete for %s", tumor$patient_id)
}


#' Preprocess a list of CamSample objects for ASM analysis
#' @param sample_list. List of CamSample objects.
#' @param config. CamConfig object.
#' @export
preprocess_asm <- function(sample_list, config) {
  for (s in sample_list) {
    # Go to next part of loop if its null
    if (is.null(s)) {
      next
    }

    # Count SNP and CpG alleles if a BAM file is provided
    cmain_asm_allele_counts(s, config)

    # Format methylation rates for ASM
    cmain_asm_make_methylation(s, config)
  }
}
