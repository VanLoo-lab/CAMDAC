#' CAMDAC tumor-normal pipeline
#' 
#' Run CAMDAC analysis on a bulk tumor and patient-matched tissue-matched tumor-adjacent normal sample.
#' 
#' @param tumor. Tumor sample data built with `create_camdac_sample()`.
#' @param normal. Normal sample data built with `create_camdac_sample()`.
#' @param config. Configuration built with `create_camdac_config()`.
#' @export
pipeline_tumor_normal <- function(tumor, normal, config){
  
  # Log
  loginfo("CAMDAC:::pipeline_tumor_normal start for %s", tumour$patient_id)
  
  # Count SNP and CpG alleles
  cmain_count_alleles(tumour, config)
  cmain_count_alleles(normal, config)
  
  # Call copy number abberations using battenberg
  cmain_make_snp_profiles(tumour, normal, config)
  cmain_run_battenberg(tumour, normal, config)
  
  # Deconvolve methylation rates
  cmain_make_methylation_profile(tumour, config)
  cmain_make_methylation_profile(normal, config)
  cmain_deconvolve_methylation(tumour, normal, config)
  
  # Call differential methylation
  cmain_call_dmps(tumour, normal, config)
  cmain_call_dmrs(tumour, config)

  # Log
  loginfo("CAMDAC:::pipeline_tumor_normal complete for %s", tumour$patient_id)
  
}