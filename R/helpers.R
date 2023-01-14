
# TODO: Better to accept columns to sort by alongside e.g. POS or start, end
#' @export
sort_genomic_dt <- function(dt, with_chr = F) {
  if (with_chr) {
    fact_levels <- paste0("chr", c(1:22, "X", "Y"))
  } else {
    fact_levels <- c(1:22, "X", "Y")
  }
  dt[, chrom := factor(chrom, levels = fact_levels)]
  return(dt[order(chrom, POS)])
}

load_camdac_opts_from_input <- function(sample_id, input_file, outdir, refdir) {
  data <- read.table(input_file) %>%
    dplyr::filter(V1 == sample_id) %>%
    unlist() %>%
    as.character()
  opt <- list()
  opt$patient_id <- data[1]
  opt$tumour_bam <- data[2]
  opt$normal_bam <- data[3]
  opt$sex <- data[4]
  opt$reference_dir <- refdir
  opt$outdir <- outdir

  tumour <- CamSample(
    patient_id = opt$patient_id,
    sex = opt$sex,
    sample_id = "T",
    sample_type = "tumour",
    bam_file = opt$tumour_bam
  )

  normal <- CamSample(
    patient_id = opt$patient_id,
    sex = opt$sex,
    sample_id = "N",
    sample_type = "normal",
    bam_file = opt$normal_bam
  )

  config <- CamConfig(
    camdac_refs = opt$reference_dir,
    outdir = opt$outdir,
    build = "hg38",
    bsseq = "wgbs",
    bsseq_lib = "pe",
    n_cores = 10,
    n_seg_split = 1000
  )

  return(list(tumour, normal, config))
}

setup_cna_inject_subdir <- function(tumour, normal, config, subdir_name) {
  # Name subdirectory
  wgs_outdir <- fs::path(config$outdir, subdir_name)

  # Sym-link allele counts and tsnps file
  allele_counts_t <- get_fpath(tumour, config, "counts")
  wgs_allele_counts_t <- fs::path(wgs_outdir, gsub(config$outdir, "", allele_counts_t))

  allele_counts_n <- get_fpath(normal, config, "counts")
  wgs_allele_counts_n <- fs::path(wgs_outdir, gsub(config$outdir, "", allele_counts_n))

  tsnps_f <- get_fpath(tumour, config, "tsnps")
  wgs_tsnps_f <- fs::path(wgs_outdir, gsub(config$outdir, "", tsnps_f))

  # Ensure parent directories exist fore ach new file
  create_parent <- function(x) fs::dir_create(fs::path_dir(x))
  create_calls <- sapply(
    c(wgs_allele_counts_t, wgs_allele_counts_n, wgs_tsnps_f),
    create_parent
  )

  # Create symlinks and update config to new output directory
  fs::file_copy(allele_counts_t, wgs_allele_counts_t, overwrite = T)
  fs::file_copy(allele_counts_n, wgs_allele_counts_n, overwrite = T)
  fs::file_copy(tsnps_f, wgs_tsnps_f, overwrite = T)

  # Update config outdir
  config$outdir <- wgs_outdir
  return(config)
}

#' Cache existing CAMDAC results into a sub-directory so that the current ones can be
#' overwritten by the refitting pipeline
#'    Decided this is unnecessary as the initial results were so wrong.
# setup_cache_refit <- function(tumour, normal, config, cache_path){
#
# }

#' Exported only for development
#' @keywords internal
helper_camdac_pileup <- function(bam_file, seg, loci_dt) {
  paired_end <- T
  drop_ccgg <- T
  bam_dt <- CAMDAC:::get_reads_in_segments(bam_file, seg, paired_end = paired_end, min_mapq = 0)
  if (nrow(bam_dt) == 0) {
    return(empty_count_alleles_result())
  }
  bam_dt <- CAMDAC:::format_bam_for_loci_overlap(bam_dt, paired_end = paired_end)
  bam_dt <- CAMDAC:::annotate_bam_with_loci(bam_dt, loci_dt, drop_ccgg = drop_ccgg, paired_end = paired_end)
  bam_dt <- CAMDAC:::drop_positions_outside_segments(bam_dt, seg)
  bam_dt <- CAMDAC:::fix_pe_overlap_at_loci(bam_dt)
  bam_dt <- CAMDAC:::add_loci_read_position(bam_dt)
  bam_dt <- CAMDAC:::fix_pe_strand_with_flags(bam_dt)
  bam_dt <- CAMDAC:::get_alleles_and_qual(bam_dt)
  bam_dt <- CAMDAC:::drop_pe_fields(bam_dt)
  bam_dt <- CAMDAC:::filter_clipped_dinucleotides(bam_dt)
  bam_dt <- CAMDAC:::annotate_nucleotide_counts(bam_dt)
  bam_dt$total_depth <- 1
  bam_dt <- CAMDAC:::get_snp_allele_counts(bam_dt)
  return(bam_dt)
}


#' Parse ASCAT and bb output directories to load CNA data
#'
#'  See "annotate_copy_number" func
#' A function required to load copy number for a tumour sample from camdac, either ascat or bb,
#'   result should be: chrom, start, end, nA, nB, CN (total), seg_min and seg_max.
#'   This should also include the purity and ploidy. As a separate list?
#'   note that seg_min and seg_max are actually duplicates of the start and end columns, required to
#'   keep track of the ascat segment positions after overalp
#'      WARN: This drops sex chromosome but not implimented. Also should drops CN=0 (hom del) regions
#' @export
load_cna_data <- function(tumour, config, data_type) {
  if (data_type == "ascat") {
    return(
      load_cna_data_ascat(tumour, config)
    )
  } else if (data_type == "battenberg") {
    return(
      load_cna_data_battenberg(tumour, config)
    )
  } else {
    stop("Invalid data type argument given")
  }
}


load_cna_data_ascat <- function(tumour, config) {
  ascat.output <- qs::qread(get_fpath(tumour, config, "ascat"))
  purity <- ascat.output$aberrantcellfraction
  ploidy <- ascat.output$ploidy
  fit <- ascat.output$goodnessOfFit

  seg <- ascat.output$segments_raw
  cna_clean <- data.table(
    chrom = factor(seg$chr, levels = c(1:22, "X", "Y")),
    start = as.numeric(as.character(seg$startpos)),
    end = as.numeric(as.character(seg$endpos)),
    nA = seg$nMajor, nB = seg$nMinor,
    CN = seg$nMajor + seg$nMinor,
    seg_min = as.numeric(as.character(seg$startpos)),
    seg_max = as.numeric(as.character(seg$endpos))
  )

  setkeyv(cna_clean, cols = c("chrom", "start", "end"))

  cna_clean$purity <- purity
  cna_clean$ploidy <- ploidy
  cna_clean$fit <- fit
  cna_clean$pipeline <- "ascat"
  return(cna_clean)
}

load_cna_data_ascat_wgs <- function(ascat_output_file) {
  load(ascat_output_file)
  purity <- ascat.output$aberrantcellfraction
  ploidy <- ascat.output$ploidy
  fit <- ascat.output$goodnessOfFit

  seg <- ascat.output$segments_raw
  cna_clean <- data.table(
    chrom = factor(seg$chr, levels = c(1:22, "X", "Y")),
    start = as.numeric(as.character(seg$startpos)),
    end = as.numeric(as.character(seg$endpos)),
    nA = seg$nMajor, nB = seg$nMinor,
    CN = seg$nMajor + seg$nMinor,
    seg_min = as.numeric(as.character(seg$startpos)),
    seg_max = as.numeric(as.character(seg$endpos))
  )

  setkeyv(cna_clean, cols = c("chrom", "start", "end"))

  result <- list(
    purity = purity, ploidy = ploidy, fit = fit, ascna = cna_clean
  )

  return(result)
}

load_cna_data_battenberg <- function(tumour, config, bb_raw = FALSE, bb_dir = NA) {
  # Allows us to use this helper for non-CAMDAC directories
  if (is.na(bb_dir)) {
    bb_dir <- fs::path_dir(get_fpath(tumour, config, "battenberg"))
  }

  # Load purity and ploidy
  pp_file <- fs::dir_ls(bb_dir, glob = "*purity_ploidy.txt*")
  pp <- data.table::fread(pp_file)
  fit_file <- fs::dir_ls(bb_dir, glob = "*rho_and_psi.txt")
  ff <- suppressWarnings(data.table::fread(fit_file))
  pp$fit <- round(ff[is.best == T]$distance * 100, 3)

  # Load cna file
  cna_file <- fs::dir_ls(bb_dir, glob = "*_subclones.txt")
  cna_file <- cna_file[!grepl("chrX", cna_file)]
  bb_cna <- data.table::fread(cna_file)
  bb_cna_fields <- c(
    "chr", "startpos", "endpos", "ntot", "nMaj1_A", "nMin1_A",
    "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A"
  )
  bb_cna <- bb_cna[, ..bb_cna_fields]
  bb_cna <- bb_cna[!is.na(nMaj1_A)]

  # Ensure Maj/Min2 are integer. Required when NA in field for compleetely clonal sample is read
  bb_cna$nMaj2_A <- as.integer(bb_cna$nMaj2_A)
  bb_cna$nMin2_A <- as.integer(bb_cna$nMin2_A)

  # Select major copy number from joint solutions
  bb_cna[, selector := (bb_cna$frac1_A > bb_cna$frac2_A) | is.na(bb_cna$frac2_A)]
  bb_cna[, nA := data.table::fifelse(
    selector, nMaj1_A, nMaj2_A
  )]
  bb_cna[, nB := data.table::fifelse(
    selector, nMin1_A, nMin2_A
  )]
  bb_cna[, CN := nA + nB]

  # Return if raw data required
  if (bb_raw) {
    return(bb_cna)
  }

  # Finalise data for export
  cna_clean <- bb_cna[, .(
    chrom = factor(chr, levels = c(1:22, "X", "Y")), start = startpos, end = endpos,
    nA, nB, CN, seg_min = startpos, seg_max = endpos
  )]
  setkeyv(cna_clean, cols = c("chrom", "start", "end"))

  cna_clean$purity <- pp$purity
  cna_clean$ploidy <- pp$ploidy
  cna_clean$fit <- pp$fit
  cna_clean$pipeline <- "battenberg"
  return(cna_clean)
}

# Helper to load clonal profile from battenberg output directory
load_clonal_bb <- function(bb_dir) {
  # Meta
  sample_id <- fs::path_file(bb_dir)

  # Load purity and ploidy
  pp_file <- fs::dir_ls(bb_dir, glob = "*purity_ploidy.txt*")
  pp <- data.table::fread(pp_file)
  fit_file <- fs::dir_ls(bb_dir, glob = "*rho_and_psi.txt")
  ff <- suppressWarnings(data.table::fread(fit_file))
  pp$fit <- round(ff[is.best == T]$distance * 100, 3)

  # Load cna file
  cna_file <- fs::dir_ls(bb_dir, glob = "*_subclones.txt")
  cna_file <- cna_file[!grepl("chrX", cna_file)]
  bb_cna <- data.table::fread(cna_file)
  bb_cna_fields <- c(
    "chr", "startpos", "endpos", "ntot", "nMaj1_A", "nMin1_A",
    "frac1_A", "nMaj2_A", "nMin2_A", "frac2_A"
  )
  bb_cna <- bb_cna[, ..bb_cna_fields]
  bb_cna <- bb_cna[!is.na(nMaj1_A)]

  # Ensure Maj/Min2 are integer. Required when NA in field for compleetely clonal sample is read
  bb_cna$nMaj2_A <- as.integer(bb_cna$nMaj2_A)
  bb_cna$nMin2_A <- as.integer(bb_cna$nMin2_A)

  # Select major copy number from joint solutions
  bb_cna[, selector := (bb_cna$frac1_A > bb_cna$frac2_A) | is.na(bb_cna$frac2_A)]
  bb_cna[, nA := data.table::fifelse(
    selector, nMaj1_A, nMaj2_A
  )]
  bb_cna[, nB := data.table::fifelse(
    selector, nMin1_A, nMin2_A
  )]
  bb_cna[, CN := nA + nB]

  # Finalise data for export
  cna_clean <- bb_cna[, .(
    chrom = factor(chr, levels = c(1:22, "X", "Y")), start = startpos, end = endpos,
    nA, nB, CN, seg_min = startpos, seg_max = endpos
  )]
  setkeyv(cna_clean, cols = c("chrom", "start", "end"))

  cna_clean$sample_id <- sample_id
  cna_clean$purity <- pp$purity
  cna_clean$ploidy <- pp$ploidy
  cna_clean$fit <- pp$fit

  return(cna_clean)
}

camdac_winsorize_tsnps <- function(tumour, config) {
  # Read data, winsorize BAF, delete outliers at 0/1 & write
  tsnps_output_file <- CAMDAC::get_fpath(tumour, config, "tsnps")
  tsnps <- data.table::fread(tsnps_output_file)
  tsnps_hets <- tsnps[between(BAFr_n, 0.15, 0.85), .(chrom, POS, BAFr)]
  baf_outliers <- winsorize(tsnps_hets$BAFr)$outliers

  # CAMDAC rule : We remove winsorize outliers that are within BAF 0 and 1
  # so that we aren't artificially removing 0s and 1s at imbalanced regions
  tsnps_hets <- tsnps_hets[baf_outliers][BAFr == 0 | BAFr == 1]

  setkey(tsnps, chrom, POS)
  tsnps <- tsnps[!tsnps_hets]
  tsnps <- sort_genomic_dt(tsnps)
  data.table::fwrite(tsnps, tsnps_output_file)
  return(tsnps_output_file)
}

chelper_import_pon_meth <- function(tumour, normal_id, config, pon_file) {
  # Imports PON file into the same patient folder as the tumour sample.
  normal_pon_tpid <- CamSample(
    tumour$patient_id, tumour$sex, normal_id, "normal", NA
  )
  outfile <- get_fpath(normal_pon_tpid, config, "methylation")
  outdir <- fs::path_dir(outfile)
  fs::dir_create(outdir)
  fs::file_copy(pon_file, outfile, overwrite = T)
  return(normal_pon_tpid) # Normal sample for deconvolution
}

read_segments_bed <- function(bed_file) {
  # Read segments from bed file
  segments <- data.table::fread(bed_file, header = FALSE)
  # Read as GRanges
  segments <- GRanges(segments$V1, IRanges(segments$V2, segments$V3))

  # Combine overlapping segments
  segments <- reduce(segments)

  # Convert to GRangesList
  segments <- split(segments, seq_len(length(segments)))

  return(segments)
}

fread_chrom <- function(x, ...) {
  # Read a file with a chrom column and ensure it is a factor
  # Used for counts and snp files where X and Y are missing
  x <- data.table::fread(x, ...)
  x$chrom <- as.character(x$chrom)
  return(x)
}
