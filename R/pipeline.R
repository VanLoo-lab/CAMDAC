#' CAMDAC analysis pipeline
#' 
#' @param tumor Tumor `CamSample` object for deconvultion.
#' @param germline Patient-matched normal `CamSample` object. May be NULL if `tumor` has CNA calls already.
#' @param infiltrates Normal `CamSample` as a proxy for infiltrating normal methylation.
#' @param origin Normal `CamSample` representing cell of origin for tumor-normal differential methylation.
#' @param config Configuration built with `CamConfig()`.
#' @export
pipeline <- function(tumor, germline, infiltrates, origin, config) {
  if (config$bsseq == "wgbs"){
    pipeline_wgbs(tumor, germline, infiltrates, origin, config)
  } else if (config$bsseq == "rrbs") {
    pipeline_rrbs(tumor, germline, infiltrates, origin, config)
  } else {
    stop("Unsupported bsseq type. Please use 'wgbs' or 'rrbs'.")
  }
}

#' Run CAMDAC WGBS analysis on a bulk tumor and patient-matched tissue-matched tumor-adjacent normal sample.
#'
#' @param tumor Tumor `CamSample` object for deconvultion.
#' @param germline Patient-matched normal `CamSample` object. May be NULL if `tumor` has CNA calls already.
#' @param infiltrates Normal `CamSample` as a proxy for infiltrating normal methylation.
#' @param origin Normal `CamSample` representing cell of origin for tumor-normal differential methylation.
#' @param config Configuration built with `CamConfig()`.
#' @export
pipeline_wgbs <- function(tumor, germline = NULL, infiltrates = NULL, origin = NULL, config) {
  # Log
  loginfo("CAMDAC:::pipeline start for %s", tumor$patient_id)

  # Preprocess CpG, SNP and methylation data for all samples
  preprocess_wgbs(
    list(tumor, germline, infiltrates, origin),
    config
  )

  # Combine tumor-germline SNPs and call CNAs
  cmain_bind_snps(tumor, germline, config)
  cmain_call_cna(tumor, config)

  # Run deconvolution
  cmain_deconvolve_methylation(tumor, infiltrates, config)

  # Call differential methylation
  cmain_call_dmps(tumor, origin, config)
  cmain_call_dmrs(tumor, config)

  # Log
  loginfo("CAMDAC WGBS pipeline complete for %s", tumor$patient_id)
}

#' Preprocess a list of CamSample objects for analysis
#' @param sample_list. List of CamSample objects.
#' @param config. CamConfig object.
#' @export
preprocess_wgbs <- function(sample_list, config) {
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


#' Call CAMDAC for a tumor and patient-matched normal sample
#' 
#' @param tumor Tumor `CamSample` object for deconvultion.
#' @param germline Patient-matched normal `CamSample` object. May be NULL if `tumor` has CNA calls already.
#' @param infiltrates Normal `CamSample` as a proxy for infiltrating normal methylation.
#' @param origin Normal `CamSample` representing cell of origin for tumor-normal differential methylation.
#' @param config Configuration built with `CamConfig()`.
#' @export
#'
#' @param patient_id character. Patient identifier
#' @param tumor_id character. Tumor sample identifier
#' @param normal_id character. Normal sample identifier
#' @param tumor_bam character. Full path to tumor bam file
#' @param normal_bam character. Full path to normal bam file
#' @param sex character. Patient sex: "XX" for female or "XY" for male
#' @param path character. Full path to CAMDAC output directory
#' @param pipeline_files character. Full path to parent directory containing CAMDAC pipeline_files
#' @param build character. Genome build: "hg19" or "hg38"
#' @param min_tumor integer. Minimum read filter for tumor samples
#' @param min_normal integer. Minimum read filter for normal samples
#' @param n_cores integer. Number of cores to use for parallel processing
#' @param mq integer. Minimum mapping quality filter
#' @export
pipeline_rrbs <- function(tumor, germline, infiltrates, origin, config){

  # Preprocess RRBS normal samples
  for (s in list(germline, infiltrates, origin)){

    # Go to next part of loop if its null
    if (is.null(s)) {
        next
    }

    preprocess_rrbs_normal(
      patient_id = s$patient_id , sample_id = s$id, bam_file = s$bam,
      min_tumor = 1, min_norm = config$min_normal_cov, mq = config$min_mapq,
      sex = s$sex, path = config$outdir,
      pipeline_files = config$refs, build = config$build,
      n_cores = config$n_cores
    )
  }

  # Main : Process RRBS tumour using the design from input files

  # Setup
  patient_id <- tumor$patient_id
  sample_id <- tumor$id
  bam_file <- tumor$bam
  sex <- tumor$sex
  path <- config$outdir
  pipeline_files <- config$refs
  build <- config$build
  n_cores <- config$n_cores
  min_tumor <- 1
  min_normal <- config$min_normal_cov
  mq <- config$min_mapq
  paired_end <- ifelse(config$lib=="pe", TRUE, FALSE)

  # Define expected ac file
  ac_file = file.path(
    path, patient_id, "Allelecounts", sample_id,
    paste0(patient_id, ".", sample_id, ".SNPs.CpGs.all.sorted.RData")
  )

  if (!file.exists(ac_file)) {
    loginfo("CAMDAC:::preprocess_rrbs_tumor: %s:%s", patient_id, sample_id)
    # Run allele counter for tumor sample
    for (a in 1:25) {
        get_allele_counts(
            i = a, patient_id = patient_id, sample_id = sample_id,
            sex = sex, bam_file = bam_file, mq = mq,
            path = path, path_to_CAMDAC = pipeline_files,
            build = build, n_cores = n_cores, test = FALSE, paired_end=paired_end
        )
    }

    # Merge allele counts
    format_output(
        patient_id, sample_id, sex, is_normal=FALSE, path, pipeline_files, build
    )

  } else {
    loginfo("CAMDAC:::preprocess_rrbs_tumor: %s already exists, skipping counts.", ac_file)
  }

  # Create SNP files and run ASCAT (tumor)
  cna_file = file.path(
        path, patient_id, "Copy_number", sample_id,
        paste0(patient_id, ".", sample_id, ".ascat.output.RData")
    )
  if (!file.exists(cna_file)){
    loginfo("CAMDAC:::ASCAT.m Tumor")
    run_ASCAT.m(
        patient_id, sample_id, sex,
        patient_matched_normal_id = germline$id,
        path, pipeline_files, build,
        min_normal, min_tumor,
        n_cores, reference_panel_coverage = NULL
    )
  } else {
      loginfo("CAMDAC:::pipeline:rrbs: %s already exists, skipping CNA", cna_file)
  }

  # Process methylation info for copy number profiling and plot summary.
  loginfo("CAMDAC:::run_methylation_data_processing Tumor")
  run_methylation_data_processing(
      patient_id, sample_id,
      normal_infiltrates_proxy_id = infiltrates$id,
      normal_origin_proxy_id = origin$id,
      path, min_normal, min_tumor, n_cores,
      reference_panel_normal_infiltrates = NULL,
      reference_panel_normal_origin = NULL
  )

  # Get purified methylation rate
  loginfo("CAMDAC:::get_pure_tumour_methylation Tumor")
  get_pure_tumour_methylation(
      patient_id = patient_id, sample_id = sample_id, sex = sex,
      normal_infiltrates_proxy_id = infiltrates$id,
      path, pipeline_files, build,
      n_cores, reseg = FALSE
  )

  # Get DMP and DMR calls
  loginfo("CAMDAC:::get_differential_methylation Tumor")
  get_differential_methylation(
      patient_id = patient_id, sample_id = sample_id, sex = sex,
      normal_origin_proxy_id = origin$id,
      path, pipeline_files, build,
      effect_size = 0.2, prob = 0.99,
      min_DMP_counts_in_DMR = 5, min_consec_DMP_in_DMR = 4,
      n_cores, reseg = FALSE, bulk = FALSE
  )

  loginfo("CAMDAC RRBS pipeline complete for %s", tumor$patient_id)
}

preprocess_rrbs_normal <- function(patient_id, sample_id, bam_file, min_tumor,
                              min_normal, mq, sex, path, pipeline_files, build, n_cores) {

    # For normals, CAMDAC-RRBS expects same ID
    normal_id = sample_id

    # Define expected allele counts.
    ac_file = file.path(
      path, patient_id, "Allelecounts", sample_id,
      paste0(patient_id, ".", sample_id, ".SNPs.CpGs.all.sorted.RData")
    )

    loginfo("CAMDAC:::preprocess_rrbs_normal: %s:%s", patient_id, sample_id)
    if(!file.exists(ac_file)) {

      # Run allele counter for normal sample
      for (a in 1:25) {
          get_allele_counts(
              i = a, patient_id = patient_id, sample_id = sample_id, sex = sex, bam_file, mq = mq,
              path, pipeline_files, build, n_cores, test = FALSE
          )
      }

      # Merge allele counts
      is_normal <- ifelse(sample_id == normal_id, TRUE, FALSE)
      format_output(
          patient_id, sample_id, sex, is_normal, path, pipeline_files, build
      )
    } else {
        loginfo("CAMDAC:::preprocess_rrbs_normal: %s already exists, skipping counts.", ac_file)
    }

    # Create SNP files (normal) or run ASCAT (tumor)
    snp_file = file.path(
        path, patient_id, "Copy_number", sample_id,
        paste0(patient_id, ".", sample_id, ".SNPs.RData")
    )
    if (!file.exists(snp_file)){
      run_ASCAT.m(
          patient_id = patient_id, sample_id = sample_id, sex = sex,
          patient_matched_normal_id = normal_id,
          path = path, path_to_CAMDAC = pipeline_files, build = build,
          min_normal = min_normal, min_tumour = NULL,
          n_cores = n_cores, reference_panel_coverage = NULL
      )
    } else {
        loginfo("CAMDAC:::preprocess_rrbs_normal: %s already exists, skipping SNP prep.", snp_file)
    }

    # Process methylation info for copy number profiling and plot summary.
    meth_file = file.path(
        path, patient_id, "Methylation", sample_id, "dt_normal_m.RData"
    )
    if (!file.exists(meth_file)){
        run_methylation_data_processing(
            patient_id, sample_id,
            normal_infiltrates_proxy_id = normal_id,
            normal_origin_proxy_id = normal_id,
            path, min_normal, min_tumor, n_cores,
            reference_panel_normal_infiltrates = NULL,
            reference_panel_normal_origin = NULL
        )
    } else {
        loginfo("CAMDAC:::preprocess_rrbs_normal: %s already exists, skipping methylation prep.", meth_file)
    }
}


preprocess_rrbs_tumor <- function(patient_id, sample_id, normal_id, bam_file, min_tumor,
                              min_normal, mq, sex, path, pipeline_files, build, n_cores) {
    # Run allele counter for normal sample
    for (a in 1:25) {
        get_allele_counts(
            i = a, patient_id = patient_id, sample_id = sample_id, sex = sex, bam_file, mq = mq,
            path, pipeline_files, build, n_cores, test = FALSE
        )
    }

    # Merge allele counts
    # Set normal status based on whether sample and normal ID match
    is_normal <- ifelse(sample_id == normal_id, TRUE, FALSE)
    format_output(
        patient_id, sample_id, sex, is_normal, path, pipeline_files, build
    )

    # Create SNP files (normal) or run ASCAT (tumor)
    run_ASCAT.m(
        patient_id, sample_id, sex,
        patient_matched_normal_id = normal_id,
        path, pipeline_files, build,
        min_normal, min_tumor,
        n_cores, reference_panel_coverage = NULL
    )

    # Process methylation info for copy number profiling and plot summary.
    run_methylation_data_processing(
        patient_id, sample_id,
        normal_infiltrates_proxy_id = normal_id,
        normal_origin_proxy_id = normal_id,
        path, min_normal, min_tumor, n_cores,
        reference_panel_normal_infiltrates = NULL,
        reference_panel_normal_origin = NULL
    )
}
