
#' Build CAMDAC sample object
#' 
#' @param patient_id An identifier for the patient
#' @param patient_sex The sex of the patient, "XX" or "XY". Required for CNA calling.
#' @param sample_type Currently unused field.
#' @param bam_file Sample BAM file
#' @export
create_camdac_sample <- function(
  patient_id, patient_sex, sample_id, sample_type, bam_file
){
  return(
    list(
      patient_id=patient_id,
      patient_sex=patient_sex,
      sample_id=sample_id,
      sample_type=sample_type, #TODO: define for tumour and multiple matched normals
      bam_file=bam_file
    )
  )
}

#' Set CAMDAC configuration
#' 
#' @param outdir A path to save CAMDAC results. The results folder structure follows the format 
#'   PATIENT/DATASET/SAMPLE/.
#'  @param bsseq Bisulfite sequencing platform. Choose between "wgbs" or "rrbs".
#'  @param bsseq_lib Bisulfite sequencing library. Choose "pe" for paired end or "se" for single end.
#'  @param build Reference genome build. Choose "hg38" or "hg19".
#'  @param n_cores Number of cores to process CAMDAC data in parallel wherever possible.
#'  @param camdac_refs Path to CAMDAC reference files. If this is not given, CAMDAC searched the
#'    environment variable CAMDAC_PIPELINE_FILES. If this is not set, CAMDAC looks in the current 
#'    working directory.
#'  @param n_seg_split In WGBS mode, the reference genome is split into segments for parallele
#'    processing. This parameter selects the number of segments per chromosome.
#'  @param min_mapq Minimum mapping quality filter for BAM file reads. Used in `cmain_allele_counts()`.
#'  @param min_cov Minimum coverage filter for: DNA methylation, Normal SNP selection.
#'  @param ascat_rho_manual ASCAT rho value for refitting
#'  @param ascat_psi_manual ASCAT/Battenberg psi value for refitting
#'  @param beaglejar Path to beagle Java file for Battenberg haplotyping (WGBS only)
#'  @export
create_camdac_config <- function(
  outdir, bsseq, bsseq_lib, build, n_cores=1, camdac_refs=NULL, n_seg_split=50,
  min_mapq=1, min_cov=3, ascat_rho_manual=NULL, ascat_psi_manual=NULL, beaglejar=NULL
){
  # TODO: Validate inputs
  
  # Create output directory if it doesn't exist and set to absolute path
  fs::dir_create(outdir)
  outdir <- fs::path_real(outdir)
  
  # Set camdac references if not they do not exist
  refs = ifelse(is.null(camdac_refs), pipeline_files(), camdac_refs)
  
  # Set beagle jar if it does no exist
  bjar = ifelse(is.null(beaglejar),
                get_reference_files(
                  list(
                    camdac_refs = refs,
                    build = build,
                    bsseq = bsseq
                  ),
                  "beagle_jar"
                ),
                beaglejar)
  
  return(list(
      camdac_refs=refs,
      outdir=outdir,
      build=build,
      bsseq=bsseq,
      bsseq_lib=bsseq_lib,
      n_cores=n_cores,
      n_seg_split=n_seg_split,
      min_mapq=min_mapq,
      min_cov=min_cov,
      ascat_rho_manual=ascat_rho_manual,
      ascat_psi_manual=ascat_psi_manual,
      beaglejar=bjar
  ))
}

is_pe <- function(config){
  # Returns TRUE if sample is paired end.
  ifelse(config$bsseq_lib=="pe", TRUE, FALSE)
}

is_ccgg <- function(config){
  # Returns TRUE if ccgg should be included in camdac run
  ifelse(config$bsseq=="wgbs", TRUE, FALSE)
}

# Create/confirm output directories

#' @export
# Get output directory for function
build_output_name <- function(sample, config, code, dir=FALSE){
  # Set output file name 
  output_name <- dplyr::case_when(
    code=="allele_counts" ~ fs::path(
      config$outdir, sample$patient_id, "Allelecounts", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "SNPs", "CpGs", "all", "sorted", "csv", "gz", sep="."
      )
    ),
    code=="segment_split" ~ fs::path( # TEST: Tempfile to place segments
      config$outdir, sample$patient_id, "Allelecounts", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "segment", "counts", fs::path_file(tempfile()), "fst", sep="."
      )
    ),
    code=="methylation" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "m", "csv", "gz", sep="."
      )
    ),
    code=="pure_methylation" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "pure", "qs", sep="."
      )
    ),
    code=="dmps" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "CAMDAC_results_per_CpG", "fst", sep="."
      )
    ),
    code=="dmrs" ~ fs::path(
      config$outdir, sample$patient_id, "Methylation", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "CAMDAC_annotated_DMRs", "fst", sep="."
      )
    ),
    code=="tsnps" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "SNPs", "csv", "gz", sep="."
      )
    ),
    code=="ascat" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$sample_id, paste(
        sample$patient_id, sample$sample_id, "ascat", "output", "qs", sep="."
      )
    ),
    code=="battenberg" ~ fs::path(
      config$outdir, sample$patient_id, "Copynumber", sample$sample_id, "battenberg", paste(
        sample$patient_id, sample$sample_id, "battenberg", "output", "qs", sep="."
      )
    )
  )
  
  if (dir){
    return(fs::path_dir(output_name))
  }
  
  return(output_name)
}

#' Split genome into segments for allele counting
#' @param segments_file An RDS file containing a GRanges object with each chromosome region to split
#' @param n_seg_split An integer to split each chromosome segment
#' @keywords internal
split_segments_gr <- function(segments_file, n_seg_split){
  segs_gr = readRDS(segments_file)
  # Split segments by n per chromosomse and combine into a single GRanges
  # Must be list and not GRangesList in order to combine downstream.
  
  #seqlengths = setNames(width(segs_gr),seqnames(segs_gr)) # Works because segs_gr is entire genome lengths
  #GenomicRanges::tileGenome(seqlengths, ntile=n_seg_split) # Can run do.call on this
  segs_tile = as.list(
    GenomicRanges::tile(segs_gr, n=n_seg_split)
  )
  segs_all = do.call(c, segs_tile)
  # Convert into list of GRanges tiles for parallelising
  segs_list = lapply(seq(length(segs_all)), function(x) segs_all[x])
  return(segs_list)  
}

#' Get CAMDAC reference files from config
#' @export
get_reference_files <- function(config, type_folder, glob=NULL){
  
  stopifnot(
    type_folder %in% c(
      "annotations", "battenberg", "gc_per_window", "loci_files", "repli_timing", "segments_files",
      "beagle_jar"
      )
  )
  
  # Select parent directory based on bsseq and build
  root = fs::dir_ls(config$camdac_refs, recurse=T, type="directory", regexp=paste0(
    config$bsseq, ".{1}", # File separator
    config$build, "$"
  ))
  
  # List requested reference files:
  refs = fs::dir_ls(fs::path(root, type_folder), glob=glob)
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
download_pipeline_files <- function(bsseq, directory=NULL){
  stopifnot(bsseq %in% c("wgbs", "rrbs", "test"))
  
  # Get download URL from CAMDAC index file
  url_index_file = system.file("extdata", "pipeline_files_urls.txt", package = "CAMDAC")
  urls = read.table(url_index_file, header=F, stringsAsFactors = F)
  names(urls) = c("bsseq", "link")
  link = urls[ urls$bsseq == bsseq, ][[2]]
  
  # Get download location
  #   If a directory is passed to the function, install there.
  if(!is.null(directory)){
    location=fs::path_expand(directory)
    fs::dir_create(location) # Ensure location exists
  } else {
  #   Else, get pipeline files location from environment variable
  #   The currect directory is used if environment variable is empty
    cpf_env = Sys.getenv("CAMDAC_PIPELINE_FILES")
    location = ifelse(cpf_env=="", ".", cpf_env)
    location = fs::path_expand(location)
  }
  
  # Ensure download directory path exists
  if(!fs::dir_exists(location)){
    fs::dir_create(location)
  }
  
  # Download pipeline files and unzip
  tf = tempfile()
  download.file(link, destfile=tf, method="wget")
  untar(tf, exdir=location)
  
  loginfo("Tempfile unpacked and deleted %s", tf)
  fs::file_delete(tf)
  
  loginfo("Pipeline files for %s downloaded to %s", bsseq, location)

}

# Return data table with CpG/SNP loci for segment
#' @export
load_loci_for_segment <- function(seg, loci_files){
  # Load CAMDAC loci object for all regions in seg
  
  # Select loci files that have chromosomes in seg
  chrom <- gsub("chr", "", unique(seqnames(seg)))
  chrom_as_loci_file_number <- dplyr::case_when(chrom=="X" ~ "23", chrom=="Y" ~ "24", TRUE ~ chrom)
  loci_file_regex = paste0(".*\\.",chrom_as_loci_file_number, ".RData", collapse="|") 
  loci_filenames = loci_files[grepl(loci_file_regex, loci_files)]
  
  # Early exit if no loci files matched
  if(length(loci_filenames)==0){return(NA)}
  
  # Load loci files as a single object
  loci_dt = data.table()
  for ( infile in loci_files ){
    load(infile) # Brings loci_subset into local environment
    loci <- data.frame(loci_subset[
      subjectHits(findOverlaps(seg, loci_subset))
      ])
    loci_dt = rbind(loci_dt, loci)
  }
  
  # Early exit if no loci are present in segment
  if(nrow(loci_dt)==0){return(NA)}
  
  # Format first column. Only works if data exists.
  names(loci_dt)[1] = "chrom"
  
  return(loci_dt)
}

pipeline_files = function(){
  pf = Sys.getenv("CAMDAC_PIPELINE_FILES")
  ifelse(pf=="", ".", pf)
}
