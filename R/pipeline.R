#' CAMDAC analysis pipeline
#' 
#' @param tumor Tumor `CamSample()` object for deconvultion.
#' @param germline Patient-matched normal `CamSample()` object. May be NULL if `tumor` has CNA calls already.
#' @param infiltrates Normal `CamSample()` as a proxy for infiltrating normal methylation.
#' @param origin Normal `CamSample()` representing cell of origin for tumor-normal differential methylation.
#' @param config Configuration built with `CamConfig()`.
#' @export
pipeline <- function(tumor, germline, infiltrates, origin, config) {
  if (config$bsseq == "wgbs"){
    logging::loginfo("WGBS analysis pipline selected.", logger="CAMDAC")
    pipeline_wgbs(tumor, germline, infiltrates, origin, config)
  } else if (config$bsseq == "rrbs") {
    logging::loginfo("RRBS analysis pipline selected.", logger="CAMDAC")
    pipeline_rrbs(tumor, germline, infiltrates, origin, config)
  } else {
    stop("Unsupported bsseq type. Please use 'wgbs' or 'rrbs'.")
  }
}

#' Run CAMDAC WGBS analysis on a bulk tumor and patient-matched tissue-matched tumor-adjacent normal sample.
#' @param tumor Tumor `CamSample` object for deconvultion.
#' @param germline Patient-matched normal `CamSample` object. May be NULL if `tumor` has CNA calls already.
#' @param infiltrates Normal `CamSample` as a proxy for infiltrating normal methylation.
#' @param origin Normal `CamSample` representing cell of origin for tumor-normal differential methylation.
#' @param config Configuration built with `CamConfig()`.
#' @keywords internal
pipeline_wgbs <- function(tumor, germline = NULL, infiltrates = NULL, origin = NULL, config) {
  # Log
  logging::loginfo("Pipeline start for %s", tumor$patient_id, logger="CAMDAC")

  # Preprocess CpG, SNP and methylation data for all samples
  sample_list <- list(
    tumor = tumor,
    germline = germline,
    infiltrates = if (identical(infiltrates, germline)) NULL else infiltrates,
    origin = if (identical(origin, germline)) NULL else origin
  )
  
  preprocess_wgbs(sample_list, config)

  # Combine tumor-germline SNPs and call CNAs
  cmain_bind_snps(tumor, germline, config)
  cmain_call_cna(tumor, config)

  # Run deconvolution
  cmain_deconvolve_methylation(tumor, infiltrates, config)

  # Call differential methylation
  cmain_call_dmps(tumor, origin, config)
  cmain_call_dmrs(tumor, config)

  # Log
  logging::loginfo("CAMDAC WGBS pipeline complete for %s", tumor$patient_id, logger="CAMDAC")
}

#' Preprocess a list of CamSample objects for analysis
#' @param sample_list. List of CamSample objects.
#' @param config. CamConfig object.
#' @export
#' @keywords internal

preprocess_wgbs <- function(sample_list, config) {
  for (s in sample_list) {
    # Go to next part of loop if its null
    if (is.null(s)) {
      next
    }

    logging::loginfo("Spawning isolated process for sample: %s", s$id, logger="CAMDAC")
    # callr::r() creates a fresh background R session for this specific block of code
    callr::r(
      func = function(sample, config) {
        library(CAMDAC)

        # 1. Setup the cluster inside the isolated process
        cl <- parallel::makeCluster(config$n_cores, type = "FORK")
        doParallel::registerDoParallel(cl)
        
        # 2. Run pipeline functions
        CAMDAC:::cmain_count_alleles(sample, config)
        CAMDAC:::cmain_make_snps(sample, config)
        CAMDAC:::cmain_make_methylation_profile(sample, config)
        
        # 3. Safely stop the cluster
        parallel::stopCluster(cl)
      }, 
      args = list(sample = s, config = config),
      show = TRUE
    )
  
  logging::loginfo("Finished sample %s. OS has reclaimed all memory.", s$id, logger="CAMDAC")
}
}


#' Call CAMDAC for a tumor and patient-matched normal sample
#' @param tumor Tumor `CamSample` object for deconvultion.
#' @param germline Patient-matched normal `CamSample` object. May be NULL if `tumor` has CNA calls already.
#' @param infiltrates Normal `CamSample` as a proxy for infiltrating normal methylation.
#' @param origin Normal `CamSample` representing cell of origin for tumor-normal differential methylation.
#' @param config Configuration built with `CamConfig()`.
#' @keywords internal
pipeline_rrbs <- function(tumor, germline, infiltrates, origin, config){

  # Preprocess RRBS normal samples
  for (s in list(germline, infiltrates, origin)){

    # Go to next part of loop if its null
    if (is.null(s)) {
        next
    }

    logging::loginfo("Preprocessing sample %s:%s", s$patient_id, s$id, logger="CAMDAC")
    preprocess_rrbs_normal(
      patient_id = s$patient_id , sample_id = s$id, bam_file = s$bam,
      min_tumor = 1, min_normal = config$min_normal_cov, mq = config$min_mapq,
      sex = s$sex, path = config$outdir,
      pipeline_files = config$refs, build = config$build,
      n_cores = config$n_cores, paired_end = is_pe(config), segments_bed=config$regions
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
  min_tumor <- config$min_cov
  min_normal <- config$min_normal_cov
  mq <- config$min_mapq
  paired_end <- is_pe(config)
  segments_bed <- config$regions

  # Define expected ac file
  ac_file = file.path(
    path, patient_id, "Allelecounts", sample_id,
    paste0(patient_id, ".", sample_id, ".SNPs.CpGs.all.sorted.RData")
  )

  if (!file.exists(ac_file)) {
    logging::loginfo("Preprocess tumour data: %s:%s", patient_id, sample_id, logger="CAMDAC")
    # Run allele counter for tumor sample
    for (a in 1:25) {
        get_allele_counts(
            i = a, patient_id = patient_id, sample_id = sample_id,
            sex = sex, bam_file = bam_file, mq = mq,
            path = path, path_to_CAMDAC = pipeline_files,
            build = build, n_cores = n_cores, test = FALSE, paired_end=paired_end, segments_bed=segments_bed
        )
    }

    # Merge allele counts
    format_output(
        patient_id, sample_id, sex, is_normal=FALSE, path, pipeline_files, build
    )

  } else {
    logging::loginfo("RRBS tumour allele count file already exists: %s.", ac_file, logger="CAMDAC")
  }

  # Create SNP files and run ASCAT (tumor)
  cna_file = file.path(
        path, patient_id, "Copy_number", sample_id,
        paste0(patient_id, ".", sample_id, ".ascat.output.RData")
    )
  if (!file.exists(cna_file)){
    logging::loginfo("ASCAT.m Tumor", logger="CAMDAC")
    run_ASCAT.m(
        patient_id, sample_id, sex,
        patient_matched_normal_id = germline$id,
        path, pipeline_files, build,
        min_normal, min_tumor,
        n_cores, reference_panel_coverage = NULL
    )
  } else {
      logging::loginfo("CNA file already exists: %s", cna_file, logger="CAMDAC")
  }

  # Process methylation info for copy number profiling and plot summary.
  logging::loginfo("Running DNA methylation processing for Tumour", logger="CAMDAC")
  run_methylation_data_processing(
      patient_id, sample_id,
      normal_infiltrates_proxy_id = infiltrates$id,
      normal_origin_proxy_id = origin$id,
      path, min_normal, min_tumor, n_cores,
      reference_panel_normal_infiltrates = NULL,
      reference_panel_normal_origin = NULL
  )

  # Get purified methylation rate
  logging::loginfo("Calculating pure tumour DNA methylation", logger="CAMDAC")
  get_pure_tumour_methylation(
      patient_id = patient_id, sample_id = sample_id, sex = sex,
      normal_infiltrates_proxy_id = infiltrates$id,
      path, pipeline_files, build,
      n_cores, reseg = FALSE
  )

  # Get DMP and DMR calls
  logging::loginfo("Get tumour differential methylation.", logger="CAMDAC")
  get_differential_methylation(
      patient_id = patient_id, sample_id = sample_id, sex = sex,
      normal_origin_proxy_id = origin$id,
      path, pipeline_files, build,
      effect_size = 0.2, prob = 0.99,
      min_DMP_counts_in_DMR = 5, min_consec_DMP_in_DMR = 4,
      n_cores, reseg = FALSE, bulk = FALSE
  )

  logging::loginfo("Pipeline complete for %s", tumor$patient_id, logger="CAMDAC")
}

preprocess_rrbs_normal <- function(patient_id, sample_id, bam_file, min_tumor,
                              min_normal, mq, sex, path, pipeline_files, build, n_cores, paired_end, segments_bed) {

    # For normals, CAMDAC-RRBS expects same ID
    normal_id = sample_id

    # Define expected allele counts.
    ac_file = file.path(
      path, patient_id, "Allelecounts", sample_id,
      paste0(patient_id, ".", sample_id, ".SNPs.CpGs.all.sorted.RData")
    )

    loginfo("CAMDAC:::preprocess_rrbs_normal: %s:%s", patient_id, sample_id)
    loginfo("Creating allele count files...")
    if(!file.exists(ac_file)) {

      # Run allele counter for normal sample
      for (a in 1:25) {
          get_allele_counts(
              i = a, patient_id = patient_id, sample_id = sample_id, sex = sex, bam_file, mq = mq,
              path, pipeline_files, build, n_cores, test = FALSE, paired_end=paired_end, segments_bed=segments_bed
          )
      }

      # Merge allele counts
      format_output(
          patient_id, sample_id, sex, is_normal = TRUE, path, pipeline_files, build
      )
      
      loginfo("Allele counting finished.")
    } else {
        loginfo("CAMDAC:::preprocess_rrbs_normal: %s already exists, skipping counts.", ac_file)
    }
    
    # Create SNP files (normal) or run ASCAT (tumor)
    loginfo("Creating SNP files...")
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
      
      loginfo("SNP files created.")
    } else {
        loginfo("CAMDAC:::preprocess_rrbs_normal: %s already exists, skipping SNP prep.", snp_file)
    }
    
    # Process methylation info for copy number profiling and plot summary.
    meth_file = file.path(
        path, patient_id, "Methylation", sample_id, "dt_normal_m.RData"
    )
    
    loginfo("Creating methylation files...")
    if (!file.exists(meth_file)){
        run_methylation_data_processing(
            patient_id, sample_id,
            normal_infiltrates_proxy_id = normal_id,
            normal_origin_proxy_id = normal_id,
            path, min_normal, min_tumor, n_cores,
            reference_panel_normal_infiltrates = NULL,
            reference_panel_normal_origin = NULL
        )
      loginfo("Mehylation files created.")
    } else {
        loginfo("CAMDAC:::preprocess_rrbs_normal: %s already exists, skipping methylation prep.", meth_file)
    }
    
}
