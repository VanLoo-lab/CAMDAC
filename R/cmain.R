#' Count alleles
#'
#' @param sample A camdac sample object
#' @param config A camac allele object
#' @export
cmain_count_alleles <- function(sample, config) {
  #  Check if outputs exist and skip if required
  output_filename <- get_fpath(sample, config, "counts")
  if (file.exists(output_filename) && !config$overwrite) {
    loginfo("Skipping allele counting for %s", paste0(sample$patient_id, ":", sample$id))
    return(output_filename)
  }

  # Load BAM regions to analyse (segments) as a list of GRanges
  if (is.null(config$regions)) {
    # Create segments across entire reference genome
    # The number of sections per chromosome is given by the config n_seg_split option.
    segments_rds <- get_reference_files(config, type = "segments_files")
    segments <- split_segments_gr(segments_rds, config$n_seg_split)
  } else {
    # Read segments BED file as a list of granges
    segments <- read_segments_bed(config$regions)
  }
  # Load SNP and CpG loci for reference genome
  loci_files <- get_reference_files(config, type = "loci_files")

  # Load sample data
  bam_file <- sample$bam
  paired_end <- is_pe(config)
  drop_ccgg <- is_ccgg(config)
  min_mapq <- config$min_mapq
  min_cov <- config$min_cov

  # Initialise parallel workers.
  doParallel::registerDoParallel(cores = config$n_cores)

  loginfo("Counting alleles for %s", paste0(sample$patient_id, ":", sample$id))
  # For each segment, load the appropriate SNP/CpG loci file segment and call allele counter in parallel
  #   Set warn=2 to ensure foreach fails if any of the parallel workers are terminated due to memory.
  #   without this option, foreach simply returns a warning and software continues
  options(warn = 2)
  tmpfiles <- foreach(seg = segments, .combine = "c") %dopar% {
    loci_dt <- load_loci_for_segment(seg, loci_files)
    ac_file <- cwrap_get_allele_counts(bam_file, seg, loci_dt, paired_end, drop_ccgg, min_mapq = min_mapq, min_cov = min_cov)
    tmp <- tempfile()
    fst::write_fst(ac_file, tmp)
    rm(loci_dt, ac_file, seg)
    gc()
    return(tmp)
  }
  options(warn = 0)

  # Combine temporary files with allele counts results into a single data table
  result <- foreach(i = tmpfiles, .combine = "rbind") %dopar% {
    fst::read_fst(i, as.data.table = T)
  }

  # Write to output file
  format_and_write_output(result, output_filename) # 2 lines, unnecessary function!
  # Delete temporary files
  foreach(i = tmpfiles) %dopar% {
    file.remove(i)
  }

  # Stop parallel workers. When running the pipeline multiple times in an R session,
  # R re-uses workers but does not clear memory. Hence large objects in foreach loops will remain.
  # TODO: When a single pipeline function has been created, then migrate doParallel calls there.
  doParallel::stopImplicitCluster()
  return(output_filename)
}



#' Make SNPs
#'
#' Format and save SNP file for CNA analysis (ASCAT or BATTENBERG)
#'
#' @param sample A camdac sample object
#' @param config A camdac config object
#' @export
cmain_make_snps <- function(sample, config) {
  output_file <- CAMDAC::get_fpath(sample, config, "snps")
  if (fs::file_exists(output_file) & !config$overwrite) {
    loginfo("Skipping SNP profile creation for %s", paste0(sample$id))
    return(output_file)
  }

  loginfo("Making SNP profile for %s", paste0(sample$id))

  # Load required reference files
  gc_refs <- get_reference_files(config, "gc_per_window")
  repli_ref <- get_reference_files(config, "repli_timing")
  loci_ref <- get_reference_files(config, "loci_files")

  # Load SNP profiles
  ac_file <- get_fpath(sample, config, "counts")
  snps <- load_snp_profile(ac_file, loci_ref)
  # Ensure SNPs sorted for ASCAT analysis
  snps <- sort_genomic_dt(snps)

  # Save tumor SNPs to output file
  fs::dir_create(fs::path_dir(output_file))
  data.table::fwrite(snps, file = output_file, compress = "gzip")

  # Return
  return(output_file)
}



#' Bind SNPs
#'
#' Combing tumor-normal SNP file for CNA analysis (ASCAT or BATTENBERG)
#'
#' @param tumor A camdac sample object
#' @param normal A camdac sample object
#' @param config A camdac config object
#' @export
cmain_bind_snps <- function(tumour, normal, config) {
  tsnps_output_file <- CAMDAC::get_fpath(tumour, config, "tsnps")
  if (fs::file_exists(tsnps_output_file) & !config$overwrite) {
    loginfo("Skipping SNP profile creation for %s", paste0(tumour$id, "&", normal$id))
    return(tsnps_output_file)
  }

  # Check previous pipeline steps have been run
  tsnps_f <- get_fpath(tumour, config, "snps")
  nsnps_f <- get_fpath(normal, config, "snps")
  if (!fs::file_exists(tsnps_f) & !fs::file_exists(nsnps_f)) {
    stop("Tumour and normal SNP profiles must be created before binding")
  }

  loginfo("Binding SNP profiles for %s", paste0(tumour$id, "&", normal$id))
  # Load required reference files
  gc_refs <- get_reference_files(config, "gc_per_window")
  repli_ref <- get_reference_files(config, "repli_timing")

  # Load SNP profiles
  tsnps <- fread_chrom(tsnps_f)
  nsnps <- fread_chrom(nsnps_f)

  # Annotate tumour SNPs
  tsnps <- annotate_normal(tsnps, nsnps, min_cov = config$min_cov)

  # Calculate LogR
  tsnps <- calculate_logr(tsnps) # Requires total_depth, total_depth_n
  # Correct LogR with GC and replication timing
  tsnps <- annotate_gc(tsnps, gc_refs, n_cores = config$n_cores) # Long-running
  tsnps <- annotate_repli(tsnps, repli_ref)
  tsnps[, LogR_corr := spline_regress_logr(LogR, GC, repli)]

  # Remove low coverage singletons (far apart from other SNPs).
  tsnps <- rm_low_cov_singletons(tsnps)

  # Ensure SNPs sorted for ASCAT analysis
  tsnps <- sort_genomic_dt(tsnps)

  # Save tumor SNPs to output file
  fs::dir_create(fs::path_dir(tsnps_output_file))
  data.table::fwrite(tsnps, file = tsnps_output_file, compress = "gzip")

  # Return
  return(tsnps_output_file)
}

#' Run ASCAT.m
#'
#' Expects SNP profiles to have been created using `cmain_make_snp_profiles`
#'
#' @param tumor A camdac sample object
#' @param normal A camdac sample object
#' @param config A camdac config object
#' @export
cmain_run_ascat <- function(tumour, normal, config) {
  # Skip if file exists and overwrite is false
  cna_output_name <- get_fpath(tumour, config, "cna")
  if (fs::file_exists(cna_output_name) & !config$overwrite) {
    loginfo("Skipping ASCAT analysis for %s", paste0(tumour$id))
    return(cna_output_name)
  }

  loginfo("Running ASCAT analysis for %s", paste0(tumour$id))

  # Setup output object and results directory
  out_obj <- get_fpath(tumour, config, "ascat")
  out_dir <- fs::dir_create(fs::path_dir(out_obj))

  # Load TSNPS
  tsnps <- fread_chrom(
    CAMDAC::get_fpath(tumour, config, "tsnps")
  )

  # Set Rho and Psi to NA if not given (required by ASCAT)
  if (!is.null(config$ascat_rho_manual) & !is.null(config$ascat_psi_manual)) {
    preset_rho <- config$ascat_rho_manual
    preset_psi <- config$ascat_psi_manual
  } else {
    preset_rho <- NA
    preset_psi <- NA
  }

  # Run ASCAT
  ascat_results <- run_ascat.m2(tumour, tsnps, outdir = out_dir, rho_manual = preset_rho, psi_manual = preset_psi)

  # Write ASCAT output files. QS used to serialise for faster read/write of WGBS data. RRBS uses .RData.
  ascat_output_name <- get_fpath(tumour, config, "ascat")
  ascat_frag_name <- gsub("output.qs", "frag.qs", ascat_output_name)
  ascat_bc_name <- gsub("output.qs", "bc.qs", ascat_output_name)

  qs::qsave(ascat_results$ascat.bc, ascat_bc_name)
  qs::qsave(ascat_results$ascat.frag, ascat_frag_name)
  qs::qsave(ascat_results$ascat.output, ascat_output_name)

  # Write CNA object to file for ease

  cna <- load_cna_data(tumor, config, "ascat")
  data.table::fwrite(cna, file = cna_output_name, sep = "\t", col.names = T, quote = F)

  return(cna_output_name)
}

#' Run battenberg
#'
#' Expects SNP profiles to have been created using `cmain_make_snp_profiles`
#'
#' @param tumor A camdac sample object
#' @param normal A camdac sample object
#' @param config A camdac config object
#' @export
cmain_run_battenberg <- function(tumour, normal, config) {
  cna_output_name <- get_fpath(tumour, config, "cna")
  if (fs::file_exists(cna_output_name) & !config$overwrite) {
    loginfo("Skipping Battenberg analysis for %s", paste0(tumour$id))
    return(cna_output_name)
  }

  # BB operates from within output directory, therefore we switch there to start and leave before ending
  currentwd <- getwd()
  outdir <- fs::dir_create(get_fpath(tumour, config, "battenberg", dir = T))
  setwd(outdir)

  # Convert CAMDAC objects to bb inputs
  tumour_prefix <- paste0(tumour$patient_id, "-", tumour$id)
  normal_prefix <- paste0(normal$patient_id, "-", normal$id)
  camdac_tumour_ac <- get_fpath(tumour, config, "counts")
  camdac_normal_ac <- get_fpath(normal, config, "counts")
  camdac_tsnps <- get_fpath(tumour, config, "tsnps")

  # Ensure camdac files exist
  stopifnot(all(sapply(
    c(camdac_tumour_ac, camdac_normal_ac, camdac_tsnps), fs::file_exists
  )))

  loginfo("Preparing WGBS allele counts for Battenberg")
  # TODO: Should we use existing SNP objects instead of ac?
  camdac_to_battenberg_allele_freqs(camdac_tumour_ac, tumour_prefix, camdac_normal_ac, normal_prefix,
    outdir,
    min_normal_depth = config$min_cov
  )

  loginfo("Preparing WGBS BAF and logR for Battenberg")
  prepare_wgbs_files <- camdac_to_battenberg_prepare_wgbs(tumour_prefix, normal_prefix, camdac_tsnps, outdir)

  loginfo("Running Battenberg for %s", paste0(tumour$id))
  # Define battenberg inputs.
  tumourname <- tumour_prefix
  normalname <- normal_prefix
  ismale <- ifelse(tumour$sex == "XY", TRUE, FALSE)

  # Setup battenberg references
  # `get_reference_files` returns files in subdirectory, so to get root we take the parent of the first file returned.
  bb_38_dir <- fs::path_dir(get_reference_files(config, "battenberg"))[[1]]
  beagleref.template <- paste0(bb_38_dir, "/beagle5/chrCHROMNAME.1kg.phase3.v5a_GRCh38nounref.vcf.gz")
  beagleplink.template <- paste0(bb_38_dir, "/beagle5/plink.chrCHROMNAME.GRCh38.map")
  problemloci <- paste0(bb_38_dir, "/probloci/probloci.txt.gz")
  imputeinfofile <- create_impute_info_file(bb_38_dir, outdir) # Created from template.

  # Set beagle software path. CAMDAC config creation fits by default.
  beaglejar <- config$beaglejar

  # Set default RHO and PSI based on config
  if (!is.null(config$ascat_rho_manual) & !is.null(config$ascat_psi_manual)) {
    use_preset_rho_psi <- T
    preset_rho <- config$ascat_rho_manual
    preset_psi <- config$ascat_psi_manual
  } else {
    use_preset_rho_psi <- F
    preset_rho <- NA
    preset_psi <- NA
  }

  # Limit number of cores to 6 to avoid battenberg memory errors.
  # TODO: Allele counts with 10 cores worked but battenberg with 10 gave OOM error. Setting nthreads to 5 for now.
  bb_cores <- ifelse(config$n_cores <= 6, config$n_cores, 6)
  min_normal_depth <- config$min_cov

  # Run battenberg
  loginfo("Running Battenberg")
  #   We could add another (optional) config parameter for battenberg cores?
  battenberg_wgbs_wrapper(tumourname, normalname, imputeinfofile, problemloci, ismale, beaglejar,
    beagleref.template, beagleplink.template,
    phasing_gamma = 2, nthreads = bb_cores,
    use_preset_rho_psi = use_preset_rho_psi, preset_rho = preset_rho,
    min_normal_depth = min_normal_depth, preset_psi = preset_psi
  )

  loginfo("Saving results")
  cna <- load_cna_data(tumor, config, "battenberg")
  data.table::fwrite(cna, file = cna_output_name, sep = "\t", col.names = T, quote = F)

  setwd(currentwd) # Return to original directory
  return(cna_output_name)
}

#' Make methylation
#'
#' Pre-process methylation from allele counts for CAMDAC deconvolution
#'
#' @param sample A camdac sample object
#' @param config A camdac config object
#' @export
cmain_make_methylation_profile <- function(sample, config) {
  loginfo("Preprocessing methylation data: %s", sample$patient_id)
  allele_counts <- data.table::fread(get_fpath(sample, config, "counts"))
  methylation <- process_methylation(allele_counts, min_meth_loci_reads = config$min_cov)
  rm(allele_counts)

  loginfo("Calculating HDI")
  hdi <- calculate_counts_hdi(methylation$M, methylation$UM, n_cores = config$n_cores)
  methylation <- cbind(methylation, hdi)
  rm(hdi)

  loginfo("Saving methylation profile: %s %s", sample$patient_id, sample$id)
  output_file <- get_fpath(sample, config, "methylation")
  fs::dir_create(fs::path_dir(output_file))
  data.table::fwrite(methylation, file = output_file)
  return(output_file)
}

#' Deconvolve methylation
#'
#' @param tumor A camdac sample object
#' @param normal A camdac sample object
#' @param config A camdac config object
#' @export
cmain_deconvolve_methylation <- function(tumour, normal, config) {
  loginfo("Combining tumour-normal methylation: %s", tumour$patient_id)
  # Load DNAme data and merge (one function)
  t_meth <- data.table::fread(get_fpath(tumour, config, "methylation"))
  n_meth <- data.table::fread(get_fpath(normal, config, "methylation"))
  meth_c <- combine_tumour_normal_methylation(t_meth, n_meth)

  loginfo("Loading CNAs: %s", tumour$patient_id)
  # Load copy number data from ascat.output and annotate CGs.
  meth_c <- annotate_cgs_with_cnas(meth_c, tumour)

  loginfo("Deconvolving DNAme: %s", tumour$patient_id)
  # Calculate m_t
  meth_c <- deconvolve_bulk_methylation(meth_c)

  # Filter: CN=0 , effective cov_t>= 3, is.na(mt-raw)
  meth_c <- filter_deconvolved_methylation(meth_c)

  loginfo("Calculating pure_tumour HDI: %s", tumour$patient_id)
  # Calculate m_t HDI # parallel - long-running function
  meth_c <- calculate_m_t_hdi(meth_c, config$n_cores)

  outfile <- get_fpath(tumour, config, "pure_methylation")
  qs::qsave(meth_c, outfile)
}

#' Call tumor-normal DMPs
#'
#' Single-sample DMP calling on CAMDAC-deconvolved data
#'
#' @param tumor A camdac sample object
#' @param normal A camdac sample object
#' @param config A camdac config object
#' @export
cmain_call_dmps <- function(tumour, normal, config) {
  loginfo("Calling DMPs")
  # Call DMPs between tumour and normal
  pmeth_file <- get_fpath(tumour, config, "pure_methylation")
  nmeth_file <- get_fpath(normal, config, "methylation")
  pmeth <- qs::qread(pmeth_file)
  nmeth <- data.table::fread(nmeth_file)

  # Ensure tumour and normal subset to the same CpGs only.
  overlaps <- findOverlaps(
    GRanges(seqnames = pmeth$chrom, ranges = IRanges(pmeth$start, pmeth$end)),
    GRanges(seqnames = nmeth$chrom, ranges = IRanges(nmeth$start, nmeth$end)),
    type = "equal"
  )
  pmeth <- pmeth[queryHits(overlaps), ]
  nmeth <- nmeth[subjectHits(overlaps), ]

  tmeth <- call_dmps(pmeth, nmeth, effect_size = 0.2, prob = 0.99, itersplit = 5e5, ncores = config$n_cores)
  tmeth_outfile <- get_fpath(tumour, config, "dmps")
  fst::write_fst(tmeth, tmeth_outfile)
}

#' Call tumor-normal DMRs
#'
#' Single-sample DMR calling on CAMDAC DMP data
#'
#' @param tumor A camdac sample object
#' @param normal A camdac sample object
#' @param config A camdac config object
#' @export
cmain_call_dmrs <- function(tumour, config) {
  loginfo("Calling DMRs")
  tmeth_outfile <- get_fpath(tumour, config, "dmps")
  tmeth_dmps <- fst::read_fst(tmeth_outfile, as.data.table = T)
  regions_file <- CAMDAC::get_reference_files(config, "annotations", "*all_regions_annotations*")
  regions_annotations <- fst::read_fst(regions_file, as.data.table = T)
  tmeth_dmrs <- call_dmrs(tmeth_dmps, regions_annotations, n_cores = config$n_cores)
  tmeth_dmrs_outfile <- get_fpath(tumour, config, "dmrs")
  fst::write_fst(tmeth_dmrs, tmeth_dmrs_outfile)
}
