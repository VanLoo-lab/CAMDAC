#' Run allele-specific methylation analysis pipeline
#' @param tumor. CamSample object for tumor sample.
#' @param germline. CamSample object for germline sample. Used for CNA calling.
#' @param infiltrates. CamSample object for infiltrating normal sample. Used for deconvolution.
#' @param origin. CamSample object for cell of origin sample. Used for differential methylation.
#' @param config. CamConfig object.
#' @export
#' @keywords internal
asm_pipeline <- function(tumor, germline = NULL, infiltrates = NULL, origin = NULL, config) {
  # Log
  loginfo("CAMDAC:::asm_pipeline start for %s", tumor$patient_id)
  sample_list <- list(tumor, germline, infiltrates, origin)

  # Checks that ASM SNPs file is available, otherwise, creates from bulk allele counts on germline
  # and attach ASM SNPs to infiltrates and origin objects if present
  cmain_asm_make_snps(tumor, germline, infiltrates, origin, config)

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
