#' Build CAMDAC sample object
#' @param id Unique identifier for the sample
#' @param sex The sex of the patient, "XX" or "XY". Required for CNA calling.
#' @param bam Sample BAM file. If not given, CAMDAC expects files linked with `attach_output`.
#' @param patient_id An identifier for the patient
#' @export
CamSample <- function(id, sex, bam = NULL, patient_id = "P") {
  return(
    list(
      id = id,
      sex = sex,
      bam = bam,
      patient_id = patient_id
    )
  )
}

#' Set CAMDAC configuration
#'
#' @param outdir A path to save CAMDAC results. The results folder structure
#' follows the format PATIENT/DATASET/SAMPLE/.
#' @param bsseq Bisulfite sequencing platform. Choose between "wgbs" or "rrbs".
#' @param lib Bisulfite sequencing library. Choose "pe" for paired end, "se" for single end.
#' @param build Reference genome build. Choose "hg38" or "hg19".
#' @param n_cores Number of cores to process CAMDAC data in parallel wherever possible.
#' @param regions A BED file with regions to restrict the analysis to
#' @param refs Path to CAMDAC reference files. If this is not given, CAMDAC searches the
#'   environment variable CAMDAC_PIPELINE_FILES. If this is not set, CAMDAC searches recursively in the current
#'   working directory.
#' @param min_mapq Minimum mapping quality filter used in `cmain_allele_counts()`.
#' @param min_cov Minimum coverage filter for: DNA methylation, Normal SNP selection.
#' @param overwrite Config to overwrite files if they already exist.
#' @param cna_caller The CNA caller to use. "ascat" or "battenberg". Default is "battenberg"
#' @param cna_settings A list of settings to pass to the CNA caller. rho, psi, java, beaglemaxmem
#' @param het_baf Numeric of length 2, the normal BAF window used to select
#'   heterozygous SNPs before CNA calling. Defaults to 0.15/0.85 for RRBS and
#'   0.08/0.92 for WGBS. Deep WGBS may warrant a tighter window, e.g. c(0.2, 0.8).
#' @export
CamConfig <- function(outdir, bsseq, lib, build, n_cores = 1, regions = NULL,
                      refs = NULL, n_seg_split = 50, min_mapq = 1, min_cov = 1, min_normal_cov=10, overwrite = FALSE,
                      cna_caller = "battenberg", cna_settings = NULL, het_baf = NULL) {
  # Create output directory if it doesn't exist and set to absolute path
  fs::dir_create(outdir)
  outdir <- fs::path_real(outdir)

  # Validate inputs
  stopifnot(cna_caller %in% c("ascat", "battenberg"))
  stopifnot(lib %in% c("pe", "se"))
  stopifnot(bsseq %in% c("wgbs", "rrbs"))

  # Resolve the heterozygous SNP BAF window for this protocol
  if (is.null(het_baf)) {
    het_baf <- HET_BAF_DEFAULTS[[bsseq]]
  }
  stopifnot(is.numeric(het_baf), length(het_baf) == 2, het_baf[1] < het_baf[2])

  # Set camdac references if not they do not exist
  refs <- ifelse(is.null(refs), pipeline_files(), fs::path_real(refs))

  # If using battenberg, validate that java is available and set beagle jar
  if (cna_caller == "battenberg") {
    check_java()
    bjar <- get_reference_files(
      list(
        refs = refs,
        build = build,
        bsseq = bsseq
      ),
      "beagle_jar"
    )
  } else {
    bjar <- NULL
  }

  # If using rrbs, CNA caller must be ASCAT
  if (bsseq == "rrbs") {
    if (cna_caller != "ascat"){
      logwarn("CNA caller must be `ascat` for RRBS data. This has now been set for the analysis.")
      cna_caller <- "ascat"
    }
  }

  if (bsseq == "wgbs"){
    if (lib != "pe"){
      logwarn(
        "WGBS data analysis is only available for paired-end samples. Stopping."
      )
      stop()
    }
  }

  return(list(
    refs = refs,
    outdir = outdir,
    build = build,
    bsseq = bsseq,
    bsseq_lib = lib,
    n_cores = n_cores,
    n_seg_split = n_seg_split,
    min_mapq = min_mapq,
    min_cov = min_cov,
    min_normal_cov = min_normal_cov,
    het_baf = het_baf,
    overwrite = overwrite,
    beaglejar = bjar,
    regions = regions,
    cna_caller = cna_caller,
    cna_settings = cna_settings
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

# Default normal BAF windows for heterozygous SNP selection, by protocol.
# WGBS is costly and so is typically sequenced at lower depth than RRBS; a wider
# window recovers more het SNPs. Deep WGBS can afford a tighter one, set via
# the `het_baf` argument to `CamConfig()`.
HET_BAF_DEFAULTS <- list(wgbs = c(0.08, 0.92), rrbs = c(0.15, 0.85))

# Resolve the het SNP BAF window, falling back to protocol defaults for config
# objects created before `het_baf` existed.
get_het_baf <- function(config) {
  if (!is.null(config$het_baf)) {
    return(config$het_baf)
  }
  HET_BAF_DEFAULTS[[config$bsseq]]
}

FPATH_CODES <- c(
  "counts", "meth", "pure", "dmps", "dmrs", "segment_split", "snps",
  "ascat", "battenberg", "tsnps", "cna", "asm_snps", "asm_counts", "asm_hap_stats",
  "asm_phase_map", "asm_meth", "asm_cna", "asm_meth_cna", "asm_meth_pure", "asm_ss_dmp", "asm_dmp"
)

# Create/confirm output directories
#' @export
get_fpath <- function(sample, config, code, dir = FALSE) {
  stopifnot(code %in% FPATH_CODES)

  # Set output file name
  output_name <- dplyr::case_when(
    code == "counts" ~ fs::path(
      config$outdir, sample$patient_id, "Allelecounts", sample$id, paste(
        sample$patient_id, sample$id, "SNPs", "CpGs", "all", "sorted", "csv", "gz",
        sep = "."
      )
    ),
    code == "segment_split" ~ fs::path( # TEST: Tempfile to place segments
      config$outdir, sample$patient_id, "Allelecounts", sample$id, paste(
        sample$patient_id, sample$id, "segment", "counts", fs::path_file(tempfile()), "fst",
        sep = "."
      )
    ),
    code == "meth" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "m", "csv", "gz",
        sep = "."
      )
    ),
    code == "pure" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "pure", "csv", "gz",
        sep = "."
      )
    ),
    code == "dmps" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "CAMDAC_results_per_CpG", "fst",
        sep = "."
      )
    ),
    code == "dmrs" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "CAMDAC_annotated_DMRs", "fst",
        sep = "."
      )
    ),
    code == "snps" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$id, paste(
        sample$patient_id, sample$id, "SNPs", "csv", "gz",
        sep = "."
      )
    ),
    code == "tsnps" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$id, paste(
        sample$patient_id, sample$id, "tnSNP", "csv", "gz",
        sep = "."
      )
    ),
    code == "ascat" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$id, "ascat", paste(
        sample$patient_id, sample$id, "ascat", "output", "qs",
        sep = "."
      )
    ),
    code == "battenberg" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$id, "battenberg", paste(
        sample$patient_id, sample$id, "battenberg", "output", "qs",
        sep = "."
      )
    ),
    code == "cna" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$id, paste(
        sample$patient_id, sample$id, "cna", "txt",
        sep = "."
      )
    ),
    code == "asm_counts" ~ fs::path(
      config$outdir, sample$patient_id, "AlleleSpecific", sample$id, paste(
        sample$patient_id, sample$id, "asm_counts", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_hap_stats" ~ fs::path(
      config$outdir, sample$patient_id, "AlleleSpecific", sample$id, paste(
        sample$patient_id, sample$id, "asm_hap_stats", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_phase_map" ~ fs::path(
      config$outdir, sample$patient_id, "AlleleSpecific", sample$id, paste(
        sample$patient_id, sample$id, "asm_phase_map", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_counts" ~ fs::path(
      config$outdir, sample$patient_id, "AlleleSpecific", sample$id, paste(
        sample$patient_id, sample$id, "asm_counts", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_snps" ~ fs::path(
      config$outdir, sample$patient_id, "AlleleSpecific", sample$id, paste(
        sample$patient_id, sample$id, "asm_snps", "txt",
        sep = "."
      )
    ),
    code == "asm_cna" ~ fs::path(
      config$outdir, sample$patient_id, "AlleleSpecific", sample$id, paste(
        sample$patient_id, sample$id, "asm_cna", "txt",
        sep = "."
      )
    ),
    code == "asm_meth" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "asm_meth", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_meth_cna" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "asm_meth_cna", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_meth_pure" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "asm_meth_pure", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_ss_dmp" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "asm_ss_dmp", "csv", "gz",
        sep = "."
      )
    ),
    code == "asm_dmp" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$id, paste(
        sample$patient_id, sample$id, "asm_dmp", "csv", "gz",
        sep = "."
      )
    )
  )

  if (dir) {
    return(fs::path_dir(output_name))
  }

  # When outpath does not exist, an uninitialised character object is given. Replace with empty string.
  if (length(output_name) == 0) {
    output_name <- ""
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
  root <- fs::dir_ls(config$refs, recurse = T, type = "directory", regexp = paste0(
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
  loginfo("Downloading pipeline files for %s analysis", bsseq)

  # Get download URL from CAMDAC index file
  url_index_file <- system.file("extdata", "pipeline_files_urls.txt", package = "CAMDAC")
  urls <- read.table(url_index_file, header = F, stringsAsFactors = F)
  names(urls) <- c("bsseq", "link")
  link <- urls[urls$bsseq == bsseq, ][[2]]

  # Get download location
  #   If a directory is passed to the function, install there.
  if (!is.null(directory)) {
    location <- fs::path_expand(directory)
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
  tf <- tempfile(tmpdir = ".")
  download.file(link, destfile = tf, method = "wget")
  loginfo("Unpacking tempfile (tar.gz): %s", tf)
  untar(tf, exdir = location)

  fs::file_delete(tf)
  loginfo("Tempfile unpacked and deleted %s", tf)

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
  for (infile in loci_filenames) {
    load(infile) # Brings loci_subset into local environment
    ol <- findOverlaps(seg, loci_subset)
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


#' Pre-load and Split Loci by Segments
#' 
#' Grouping by chromosome ensures each .RData file is read only once, 
#' Reduce I/O contention before parallel processing.
#'
#' @param segments A list of GRanges objects.
#' @return A list of keyed data.tables matching the length and order of 'segments'.
prepare_loci_list <- function(segments, loci_files, drop_ccgg = FALSE) {
  logging::loginfo("Pre-loading and splitting loci for parallel workers...", logger="CAMDAC")
  
  # Create an empty list to hold results
  loci_list <- vector("list", length(segments))
  # A segment may hold several ranges, and those ranges may sit on more than one
  # chromosome, so index each segment under every chromosome it touches.
  seg_chroms <- lapply(segments, function(x) unique(as.character(seqnames(x))))
  all_chroms <- unique(unlist(seg_chroms))

  for (ch in all_chroms) {
    seg_idx <- which(vapply(seg_chroms, function(cs) ch %in% cs, logical(1)))
    if (length(seg_idx) == 0) {
      next
    }

    # 1. Map chromosome name to CAMDAC file numbering
    ch_num <- dplyr::case_when(
      ch == "chrX" | ch == "X" ~ "23", 
      ch == "chrY" | ch == "Y" ~ "24", 
      TRUE ~ gsub("chr", "", ch)
    )
    
    # 2. Identify and load the specific Chromosome file
    # Uses $ to ensure exact match (e.g., .1.RData not .11.RData)
    ch_file <- loci_files[grepl(paste0("\\.", ch_num, "\\.RData$"), loci_files)]
    
    if (length(ch_file) == 0) {
      logging::logwarn("No loci file found for chromosome: %s", ch, logger="CAMDAC")
      next
    }
    
    # Load 'loci_subset' into a temporary environment to avoid name collisions
    tmp_env <- new.env()
    load(ch_file, envir = tmp_env)
    
    if (!exists("loci_subset", envir = tmp_env)) {
      logging::logerror("Object 'loci_subset' not found in %s", ch_file, logger="CAMDAC")
      next
    }
    
    # 3. Convert to keyed data.table once per chromosome
    ch_loci_dt <- data.table::as.data.table(tmp_env$loci_subset)
    data.table::setnames(ch_loci_dt, 1, "chrom")
    # Ensure chrom is character and keys are set ONCE per chromosome
    ch_loci_dt[, chrom := as.character(chrom)]
    essential_cols <- c("chrom", "start", "end", "width", "strand", "POS", "ref", "alt")
    ch_loci_dt <- ch_loci_dt[, ..essential_cols]
    
    # 4. Pre-apply Filters for WGBS
    if (drop_ccgg) {
      ch_loci_dt <- ch_loci_dt[width != 4]
    }

    # 5. SET KEY ONCE (Crucial for inherited keys in workers)
    data.table::setkey(ch_loci_dt, chrom, start, end)

    # 6. Extract loci for each segment on this chromosome.
    # Select by genomic overlap, as load_loci_for_segment() did. A containment
    # test would silently drop loci straddling a segment boundary, and would
    # recycle rather than error on segments holding more than one range.
    ch_loci_gr <- GRanges(
      seqnames = ch_loci_dt$chrom,
      ranges = IRanges(ch_loci_dt$start, ch_loci_dt$end)
    )

    for (idx in seg_idx) {
      s <- segments[[idx]]
      # Restrict to the ranges of this segment that lie on the current chromosome
      s_ch <- s[as.character(seqnames(s)) == ch]
      hits <- sort(unique(subjectHits(findOverlaps(s_ch, ch_loci_gr))))
      if (length(hits) == 0) {
        next
      }
      part <- ch_loci_dt[hits]
      # A segment spanning chromosomes accumulates loci across iterations
      loci_list[[idx]] <- if (is.null(loci_list[[idx]])) {
        part
      } else {
        rbind(loci_list[[idx]], part)
      }
    }

    # Cleanup chromosome-level data to keep Master Process lean
    rm(tmp_env, ch_loci_dt, ch_loci_gr); gc()
  }

  # Mark empty segments and key the rest, ready for zero-copy sharing in workers
  for (i in seq_along(loci_list)) {
    if (is.null(loci_list[[i]]) || nrow(loci_list[[i]]) == 0) {
      loci_list[[i]] <- NA # Explicitly mark empty
    } else {
      data.table::setkeyv(loci_list[[i]], c("chrom", "start", "end"))
    }
  }

  return(loci_list)
}

pipeline_files <- function() {
  pf <- Sys.getenv("CAMDAC_PIPELINE_FILES")
  ifelse(pf == "", fs::path_real("."), pf)
}

check_java <- function() {
  java_found <- system2("java", "-version", stderr = F, stdout = F) == 0
  if (!java_found) {
    stop(paste0(
      "Java not found. Please install Java to use Battenberg,",
      "otherwise, run set cna_caller to ASCAT and try again."
    ))
  }
}
