
#' Build CAMDAC sample object
#'
#' @param patient_id An identifier for the patient
#' @param patient_sex The sex of the patient, "XX" or "XY". Required for CNA calling.
#' @param sample_type Currently unused field.
#' @param bam_file Sample BAM file
#' @param segments_bed Optional. a BED file of genomic regions to analyse with CAMDAC.  If empty,
#'   the whole genome is processed.
#' @param cna Optional. Pre-determined allele-specific copy number data for the tumor sample. Format
#'   is a list with three entries:
#'  - `purity` fraction denoting tumor purity
#'  - `ploidy` overal tumor ploidy e.g. 3.25
#'  - `ascna` a data frame where each row is a copy number segment. Data is taken by column position
#'     where 1=chromosome, 2=start, 3=end, 4=major copy number and 5=minor copy number. If the ascna
#'     entry is left blank, the purity and ploidy entries will be used for re-fitting
#' @export
CamSample <- function(patient_id, patient_sex, sample_id, sample_type, bam_file, segments_bed = NULL, cna = NULL) {
  return(
    list(
      patient_id = patient_id,
      patient_sex = patient_sex,
      sample_id = sample_id,
      sample_type = sample_type, # TODO: define for tumour and multiple matched normals
      bam_file = bam_file,
      segments_bed = segments_bed,
      cna = validate_cna(cna)
    )
  )
}

#' Set CAMDAC configuration
#'
#' @param outdir A path to save CAMDAC results. The results folder structure follows the format
#'   PATIENT/DATASET/SAMPLE/.
#'  @param bsseq Bisulfite sequencing platform. Choose between "wgbs" or "rrbs".
#'  @param lib Bisulfite sequencing library. Choose "pe" for paired end, "se" for single end.
#'  @param build Reference genome build. Choose "hg38" or "hg19".
#'  @param n_cores Number of cores to process CAMDAC data in parallel wherever possible.
#'  @param camdac_refs Path to CAMDAC reference files. If this is not given, CAMDAC searched the
#'    environment variable CAMDAC_PIPELINE_FILES. If this is not set, CAMDAC looks in the current
#'    working directory.
#'  @param min_mapq Minimum mapping quality filter used in `cmain_allele_counts()`.
#'  @param min_cov Minimum coverage filter for: DNA methylation, Normal SNP selection.
#'  @param overwrite Config to overwrite files if they already exist.
#'  @export
CamConfig <- function(outdir, bsseq, lib, build, n_cores = 1,
                      camdac_refs = NULL, min_mapq = 1, min_cov = 3, overwrite = FALSE) {
  # TODO: Validate inputs
  # TODO: Ensure overwrite is used by all cmain* pipeline functions.

  # Create output directory if it doesn't exist and set to absolute path
  fs::dir_create(outdir)
  outdir <- fs::path_real(outdir)

  # Set camdac references if not they do not exist
  refs <- ifelse(is.null(camdac_refs), pipeline_files(), camdac_refs)

  # Set beagle jar if it does no exist
  bjar <- get_reference_files(
    list(
      camdac_refs = refs,
      build = build,
      bsseq = bsseq
    ),
    "beagle_jar"
  )

  return(list(
    camdac_refs = refs,
    outdir = outdir,
    build = build,
    bsseq = bsseq,
    bsseq_lib = lib,
    n_cores = n_cores,
    n_seg_split = 50,
    min_mapq = min_mapq,
    min_cov = min_cov,
    overwrite = overwrite,
    beaglejar = bjar
  ))
}

is_pe <- function(config) {
  # Returns TRUE if sample is paired end.
  ifelse(config$bsseq_lib == "pe", TRUE, FALSE)
}

is_ccgg <- function(config) {
  # Returns TRUE if ccgg should be included in camdac run
  ifelse(config$bsseq == "wgbs", TRUE, FALSE)
}

# Create/confirm output directories

#' @export
# Get output directory for function
build_output_name <- function(sample, config, code, dir = FALSE) {
  # Set output file name
  output_name <- dplyr::case_when(
    code == "allele_counts" ~ fs::path(
      config$outdir, sample$patient_id, "Allelecounts", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "SNPs", "CpGs", "all", "sorted", "csv", "gz",
        sep = "."
      )
    ),
    code == "segment_split" ~ fs::path( # TEST: Tempfile to place segments
      config$outdir, sample$patient_id, "Allelecounts", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "segment", "counts", fs::path_file(tempfile()), "fst",
        sep = "."
      )
    ),
    code == "methylation" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "m", "csv", "gz",
        sep = "."
      )
    ),
    code == "pure_methylation" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "pure", "qs",
        sep = "."
      )
    ),
    code == "dmps" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "CAMDAC_results_per_CpG", "fst",
        sep = "."
      )
    ),
    code == "dmrs" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "CAMDAC_annotated_DMRs", "fst",
        sep = "."
      )
    ),
    code == "tsnps" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "SNPs", "csv", "gz",
        sep = "."
      )
    ),
    code == "ascat" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "ascat", "output", "qs",
        sep = "."
      )
    ),
    code == "battenberg" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$sample_id, "battenberg", paste(
        sample$patient_id, sample$sample_id, "battenberg", "output", "qs",
        sep = "."
      )
    )
  )

  if (dir) {
    return(fs::path_dir(output_name))
  }

  return(output_name)
}

#' Split genome into segments for allele counting
#' @param segments_file An RDS file containing a GRanges object with each chromosome region to split
#' @param n_seg_split An integer to split each chromosome segment
#' @keywords internal
split_segments_gr <- function(segments_file, n_seg_split) {
  segs_gr <- readRDS(segments_file)
  # Split segments by n per chromosomse and combine into a single GRanges
  # Must be list and not GRangesList in order to combine downstream.

  # seqlengths = setNames(width(segs_gr),seqnames(segs_gr)) # Works because segs_gr is entire genome lengths
  # GenomicRanges::tileGenome(seqlengths, ntile=n_seg_split) # Can run do.call on this
  segs_tile <- as.list(
    GenomicRanges::tile(segs_gr, n = n_seg_split)
  )
  segs_all <- do.call(c, segs_tile)
  # Convert into list of GRanges tiles for parallelising
  segs_list <- lapply(seq(length(segs_all)), function(x) segs_all[x])
  return(segs_list)
}

#' Get CAMDAC reference files from config
#' @export
get_reference_files <- function(config, type_folder, glob = NULL) {
  stopifnot(
    type_folder %in% c(
      "annotations", "battenberg", "gc_per_window", "loci_files", "repli_timing", "segments_files",
      "beagle_jar"
    )
  )

  # Select parent directory based on bsseq and build
  root <- fs::dir_ls(config$camdac_refs, recurse = T, type = "directory", regexp = paste0(
    config$bsseq, ".{1}", # File separator
    config$build, "$"
  ))

  # List requested reference files:
  refs <- fs::dir_ls(fs::path(root, type_folder), glob = glob)
  return(refs)
}

#' Download CAMDAC pipeline files
#'
#' @description CAMDAC pipeline files are required for analysis. This function downloads the files to
#' the output directory and unpacks them. By default, CAMDAC searches for the files in the
#' environment variable CAMDAC_PIPELINE_FILES. If this is missing, the current directory is used.
#' @param assay Sequencing assay. Either wgbs or rrbs.
#' @param directory Optional. Directory to download files to.
#' @export
download_pipeline_files <- function(bsseq, directory = NULL) {
  stopifnot(bsseq %in% c("wgbs", "rrbs", "test"))

  # Get download URL from CAMDAC index file
  url_index_file <- system.file("extdata", "pipeline_files_urls.txt", package = "CAMDAC")
  urls <- read.table(url_index_file, header = F, stringsAsFactors = F)
  names(urls) <- c("bsseq", "link")
  link <- urls[urls$bsseq == bsseq, ][[2]]

  # Get download location
  #   If a directory is passed to the function, install there.
  if (!is.null(directory)) {
    location <- fs::path_expand(directory)
    fs::dir_create(location) # Ensure location exists
  } else {
    #   Else, get pipeline files location from environment variable
    #   The currect directory is used if environment variable is empty
    cpf_env <- Sys.getenv("CAMDAC_PIPELINE_FILES")
    location <- ifelse(cpf_env == "", ".", cpf_env)
    location <- fs::path_expand(location)
  }

  # Ensure download directory path exists
  if (!fs::dir_exists(location)) {
    fs::dir_create(location)
  }

  # Download pipeline files and unzip
  tf <- tempfile()
  download.file(link, destfile = tf, method = "wget")
  untar(tf, exdir = location)

  loginfo("Tempfile unpacked and deleted %s", tf)
  fs::file_delete(tf)

  loginfo("Pipeline files for %s downloaded to %s", bsseq, location)
}

# Return data table with CpG/SNP loci for segment
#' @export
load_loci_for_segment <- function(seg, loci_files) {
  # Load CAMDAC loci object for all regions in seg

  # Select loci files that have chromosomes in seg
  chrom <- gsub("chr", "", unique(seqnames(seg)))
  chrom_as_loci_file_number <- dplyr::case_when(chrom == "X" ~ "23", chrom == "Y" ~ "24", TRUE ~ chrom)
  loci_file_regex <- paste0(".*\\.", chrom_as_loci_file_number, ".RData", collapse = "|")
  loci_filenames <- loci_files[grepl(loci_file_regex, loci_files)]

  # Early exit if no loci files matched
  if (length(loci_filenames) == 0) {
    return(NA)
  }

  # Load loci files as a single object
  loci_dt <- data.table()
  for (infile in loci_files) {
    load(infile) # Brings loci_subset into local environment
    loci <- data.frame(loci_subset[
      subjectHits(findOverlaps(seg, loci_subset))
    ])
    loci_dt <- rbind(loci_dt, loci)
  }

  # Early exit if no loci are present in segment
  if (nrow(loci_dt) == 0) {
    return(NA)
  }

  # Format first column. Only works if data exists.
  names(loci_dt)[1] <- "chrom"

  return(loci_dt)
}

pipeline_files <- function() {
  pf <- Sys.getenv("CAMDAC_PIPELINE_FILES")
  ifelse(pf == "", ".", pf)
}

validate_cna <- function(cna) {
  if (is.null(cna)) {
    return(NULL)
  }

  names_present <- all(c("purity", "ploidy", "ascna") %in% names(cna))
  if (!names_present) {
    logerror("CNA object must have purity, ploidy and ascna fields")
    stop()
  }
  numeric_pp <- all(is.numeric(c(cna$purity, cna$ploidy)))
  if (!numeric_pp) {
    logerror("CNA purity and ploidy must be numeric")
    stop()
  }
  if (!is.null(cna$ascna)) {
    names(cna$ascna)[1:5] <- c("chrom", "start", "end", "nA", "nB")
  }
  return(cna)
}
