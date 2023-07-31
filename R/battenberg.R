
# Convert CAMDAC allele counts object to per-chromosome allele frequency files for battenberg
create_allele_frequencies_chr <- function(tsnps, outdir, sample_af_prefix, min_depth, is_tumor=T) {
  # Load data table of SNP positions
  dt <- tsnps[, .(chrom, POS, ref, alt)]
  
  if(is_tumor){
    # For tumor, alt counts are derived from BAF
    dt$alt_counts = round(tsnps$BAF * tsnps$total_counts, 0)
    dt$ref_counts = tsnps$total_counts - dt$alt_counts
  }else{
    # For normal, alt counts are derived from BAF_n
    dt$alt_counts = round(tsnps$BAF_n * tsnps$total_counts_n, 0)
    dt$ref_counts = tsnps$total_counts_n - dt$alt_counts
  }
  

  # Filter for min depth
  dt <- dt[ref_counts + alt_counts >= min_depth]
  # Split per chromosome for parallel processing
  dt_per_chrom <- split(dt, dt$chrom)
  rm(dt)

  af_files <- foreach(dt_chrom = dt_per_chrom) %dopar% {
    
    # Set output filename with correct suffix for X and Y chromosomes
    chromosome <- unique(dt_chrom$chrom)[[1]]
    outfile_chrom <- data.table::fcase(
      chromosome == "X", "23",
      chromosome == "Y", "24",
      !(chromosome %in% c("X", "Y")), chromosome
    )
    outfile <- fs::path(outdir, sprintf("%s_alleleFrequencies_chr%s.txt", sample_af_prefix, outfile_chrom))

    # Set CHR column to allelefreq format expected
    data.table::setnames(dt_chrom, "chrom", "#CHR")

    # Set counts based on ref and alt
    nucs = c("A", "C", "G", "T")
    for( nn in nucs ){
      dt_chrom[ ref == nn, paste0("Count_", nn) := ref_counts ]
      dt_chrom[ alt == nn, paste0("Count_", nn) := alt_counts ]
    }

    # Cleanup data table
    dt_chrom = dt_chrom[!is.na(POS)]
    dt_chrom[, names(dt_chrom) := lapply(.SD, function(x) {x[is.na(x)] <- 0 ; x}) ]

    # Stack ref and alt count columns
    headings <- c("#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T")
    dt_chrom = dt_chrom[, ..headings]
    dt_chrom$Good_depth = rowSums(dt_chrom[, .(Count_A, Count_C, Count_G, Count_T)])

    # Order by position. Possible as only a single chromosome is present
    dt_chrom <- dt_chrom[order(POS)]

    # Write to output file
    write.table(dt_chrom, file = outfile, sep = "\t", row.names = F, col.names = T, quote = F)
    return(outfile)
  }
  return(af_files)
}

# AF files for tumour and normal must refer to the same SNP loci
# Note that this is true of BAF and LogR but already handled by the use of tsnps object
filter_allele_frequencies_to_overlap <- function(outdir, tumour_prefix, normal_prefix) {
  chroms <- c(1:23)
  doParallel::registerDoParallel(cores = 10)
  foreach(chrom = chroms) %dopar% {
    tumour_af <- fs::path(outdir, paste0(tumour_prefix, "_alleleFrequencies_chr", chrom, ".txt"))
    normal_af <- fs::path(outdir, paste0(normal_prefix, "_alleleFrequencies_chr", chrom, ".txt"))
    chrom_files_exist <- (all(sapply(c(tumour_af, normal_af), fs::file_exists)))

    # Skip chromosomes where data is missing from either tumour or normal
    if (!chrom_files_exist) {
      return(NULL)
    }

    tumour_dt <- fread(tumour_af)
    setkey(tumour_dt, "#CHR", "POS")
    normal_dt <- fread(normal_af)
    setkey(normal_dt, "#CHR", "POS")

    stopifnot(unique(tumour_dt$`#CHR`) == unique(normal_dt$`#CHR`)) # Ensure same chroms

    # Get overlapping positions
    positions <- intersect(tumour_dt$POS, normal_dt$POS)
    tumour_dt <- tumour_dt[POS %in% positions]
    normal_dt <- normal_dt[POS %in% positions]

    # Overwrite allele freq files
    write.table(tumour_dt, file = tumour_af, sep = "\t", row.names = F, col.names = T, quote = F)
    write.table(normal_dt, file = normal_af, sep = "\t", row.names = F, col.names = T, quote = F)
  }
  doParallel::stopImplicitCluster()
}

make_allele_frequencies_chrX_chrY <- function(outdir, tumour_prefix, normal_prefix) {
  # For some reason, BB can read imputeinfo and fail to capture chrom X and Y allelfReq files, so I'm copying here.
  mapping <- list("chr23" = "chrX", "chr24" = "chrY")
  af_files <- c(fs::dir_ls(outdir, glob = "*chr23*"), fs::dir_ls(outdir, glob = "*chr24*"))
  for (i in af_files) {
    old_file <- i
    for (n in seq_along(mapping)) { # Set based on any match in mapping
      i <- gsub(names(mapping[n]), mapping[n][[1]], i)
    }
    fs::file_copy(old_file, i, overwrite = T)
  }
}

create_impute_info_file <- function(bb_38_dir, outdir) {
  impute_template <- fs::path(bb_38_dir, "imputation", "impute_info.txt")
  data <- read.table(impute_template, header = F)
  path_sub <- function(x, prefix) gsub("\\$\\{REF_PATH\\}", fs::path(bb_38_dir, prefix), x)
  data[, 2] <- path_sub(data[[2]], "shapeit2")
  data[, 3] <- path_sub(data[[3]], "imputation")
  data[, 4] <- path_sub(data[[4]], "shapeit2")

  impute_info_out <- fs::path(outdir, "impute_info.txt")
  write.table(data, file = impute_info_out, quote = F, sep = "\t", row.names = F, col.names = F)
  return(impute_info_out)
}

camdac_to_battenberg_allele_freqs <- function(tsnps, tumour_prefix, normal_prefix, outdir, min_normal_depth) {
  create_allele_frequencies_chr(tsnps, outdir, tumour_prefix, min_depth = 1, is_tumor=T)
  create_allele_frequencies_chr(tsnps, outdir, normal_prefix, min_depth = min_normal_depth, is_tumor=F)
  filter_allele_frequencies_to_overlap(outdir, tumour_prefix, normal_prefix)
  make_allele_frequencies_chrX_chrY(outdir, tumour_prefix, normal_prefix)
}

# Battenberg LogR writer taken from : https://github.com/Wedge-lab/battenberg/blob/c257a710d88b23986f936e0b7b38131279d07f7e/R/prepare_wgs.R#L348
# Battenberg BAF writer taken from : https://github.com/Wedge-lab/battenberg/blob/c257a710d88b23986f936e0b7b38131279d07f7e/R/prepare_wgs.R#L130
create_logr_and_baf_files <- function(tsnps_file, tumour_prefix, normal_prefix, outdir) {
  # Set output file column var
  outfile_columns <- c("Chromosome", "Position")

  # Load and format TSNPs file containing CAMDAC LogR+BAF
  tsnps <- fread_chrom(tsnps_file)

  data.table::setnames(tsnps, "chrom", outfile_columns[[1]])
  data.table::setnames(tsnps, "POS", outfile_columns[[2]])

  # Create mutant LogR and BAF files
  mutant_logr_outfile <- fs::path(outdir, sprintf("%s_mutantLogR.tab", tumour_prefix))
  mutant_baf_outfile <- fs::path(outdir, sprintf("%s_mutantBAF.tab", tumour_prefix))
  outfile_columns[[3]] <- tumour_prefix # Temporarily set last outfile column to mutant name for writing
  data.table::setnames(tsnps, "LogR", outfile_columns[[3]])
  readr::write_tsv(x = tsnps[, ..outfile_columns], mutant_logr_outfile)
  data.table::setnames(tsnps, outfile_columns[[3]], "LogR") # Return tnsps table value
  data.table::setnames(tsnps, "BAFr", outfile_columns[[3]])
  write.table(tsnps[, ..outfile_columns], file = mutant_baf_outfile, row.names = F, quote = F, sep = "\t", col.names = outfile_columns)
  data.table::setnames(tsnps, outfile_columns[[3]], "BAFr") # Return tnsps table value

  # Create mutant LogR gc corrected file
  mutant_logr_gc_outfile <- fs::path(outdir, sprintf("%s_mutantLogR_gcCorrected.tab", tumour_prefix))
  outfile_columns[[3]] <- tumour_prefix # Temporarily set last outfile column to mutant name for writing
  data.table::setnames(tsnps, "LogR_corr", outfile_columns[[3]])
  readr::write_tsv(x = tsnps[, ..outfile_columns], mutant_logr_gc_outfile)
  data.table::setnames(tsnps, outfile_columns[[3]], "LogR_corr") # Return tnsps table value

  # Create normal LogR and BAF files
  normal_logr_outfile <- fs::path(outdir, sprintf("%s_normalLogR.tab", normal_prefix))
  normal_baf_outfile <- fs::path(outdir, sprintf("%s_normalBAF.tab", normal_prefix))
  outfile_columns[[3]] <- normal_prefix # Temporarily set last outfile column to normal name for writing
  tsnps[, normalLogR := 0] # Set normal LogR to 0
  data.table::setnames(tsnps, "normalLogR", outfile_columns[[3]])
  readr::write_tsv(x = tsnps[, ..outfile_columns], normal_logr_outfile)
  data.table::setnames(tsnps, outfile_columns[[3]], "normalLogR") # Return tnsps table value
  data.table::setnames(tsnps, "BAFr_n", outfile_columns[[3]])
  write.table(tsnps[, ..outfile_columns], file = normal_baf_outfile, row.names = F, quote = F, sep = "\t", col.names = outfile_columns)
  data.table::setnames(tsnps, outfile_columns[[3]], "BAFr_n") # Return tnsps table value

  return(NULL)
}

#' Generate alleleCounter file from CAMDAC
#'
#' `camdac_to_battenberg_prepare_wgbs` converts CAMDAC allele counter results to a format for processing.
#'
#' @param camdac_tumour_ac CAMDAC tumour allele counts filepath. Expected *.gz
#' @param camdac_normal_ac CAMDAC normal allele couts filepath. Expected *.gz
#' @param camdac_tnsps CAMDAC tumour-normal-snps object. Expected *.gz
#' @param output_file allelecounter formatted-file output directory.
#'
#' @returns File handle for allele counter file generated
#' @keywords internal
camdac_to_battenberg_prepare_wgbs <- function(tumour_prefix, normal_prefix, camdac_tsnps, outdir) {
  # Return record of files if all files exist
  filenames <- c(
    fs::path(outdir, paste0(tumour_prefix, "_mutantLogR.tab")),
    fs::path(outdir, paste0(tumour_prefix, "_mutantLogR_gcCorrected.tab")),
    fs::path(outdir, paste0(tumour_prefix, "_mutantBAF.tab")),
    fs::path(outdir, paste0(normal_prefix, "_normalLogR.tab")),
    fs::path(outdir, paste0(normal_prefix, "_normalBAF.tab"))
  )

  if (all(sapply(filenames, fs::file_exists))) {
    return(filenames)
  }

  # Create the mutantLogR and mutantBAF files
  # Create the normalLogR and normalBAF files
  # Create the mutantLogR_gcCorrected file
  create_logr_and_baf_files(camdac_tsnps, tumour_prefix, normal_prefix, outdir)

  # TODO: Create ASCAT plots and simplify the tables by passing sample name to the write output function:
  # see: https://github.com/Wedge-lab/battenberg/blob/43686673566cf5adbd8d00e2450d70eced27696d/R/prepare_wgs.R#L154

  return(filenames)
}

check_callChrXSubclones <- function(TUMOURNAME) {
  # Helper function to check whether we can run callChrX subclones i.e. enough SNPs in non-par regions
  # N.B. May not be valid if reference significantly changes.
  PCFinput <- data.frame(Battenberg::read_table_generic(paste0(TUMOURNAME, "_mutantLogR_gcCorrected.tab")), stringsAsFactors = F)
  PCFinput <- PCFinput[which(PCFinput$Chromosome == "X" & PCFinput$Position > 2.6e6 & PCFinput$Position < 156e6), ] # get nonPAR
  colnames(PCFinput)[3] <- TUMOURNAME
  nrow(PCFinput) > 0
}

# Run phasing to end of BB pipeline
# The Battenberg::battenberg function runs the main pipeline, however in the recent dev version,
# internally parallelised tasks fail on our system due to the use of parallell::makeCluster to
# initialise foreach.
# We redefine the Battenberg wrapper to resolve this issue.
battenberg_wgbs_wrapper <- function(tumourname,
                                    normalname,
                                    imputeinfofile,
                                    problemloci,
                                    ismale,
                                    beaglejar,
                                    beagleref.template,
                                    beagleplink.template,
                                    nthreads = 5,
                                    externalhaplotypefile = NA,
                                    allelecounts_file = NULL,
                                    sampleidx = 1,
                                    usebeagle = TRUE,
                                    beaglewindow = 40,
                                    beagleoverlap = 4,
                                    beaglemaxmem = 10,
                                    beaglenthreads = 1,
                                    gccorrectprefix = NULL,
                                    data_type = "wgs",
                                    impute_exe = NULL,
                                    platform_gamma = 1,
                                    phasing_gamma = 1,
                                    segmentation_gamma = 10,
                                    segmentation_kmin = 3,
                                    phasing_kmin = 1,
                                    clonality_dist_metric = 0,
                                    ascat_dist_metric = 1,
                                    min_ploidy = 1.6,
                                    max_ploidy = 4.8,
                                    min_rho = 0.1,
                                    min_goodness = 0.63,
                                    uninformative_BAF_threshold = 0.51,
                                    min_normal_depth = 10,
                                    min_base_qual = 20,
                                    min_map_qual = 35,
                                    calc_seg_baf_option = 3,
                                    skip_allele_counting = T, # T
                                    skip_preprocessing = T, # T
                                    skip_phasing = F,
                                    javajre = "java",
                                    write_battenberg_phasing = T,
                                    multisample_relative_weight_balanced = 0.25,
                                    multisample_maxlag = 100,
                                    segmentation_gamma_multisample = 5,
                                    snp6_reference_info_file = NA,
                                    heterozygousFilter = "none",
                                    prior_breakpoints_file = NULL,
                                    GENOMEBUILD = "hg38",
                                    use_preset_rho_psi = F, # Added to expose manual setting
                                    preset_rho = NA,
                                    preset_psi = NA) {
  # Battenberg WGBS (currently) begins from haplotyping step. First, ensure expected files are present.
  stopifnot(all(sapply(
    c(
      imputeinfofile,
      beaglejar,
      problemloci,
      paste0(tumourname, c("_mutantBAF.tab", "_mutantLogR.tab")),
      paste0(normalname, c("_normalBAF.tab", "_normalLogR.tab"))
    ),
    fs::file_exists
  )))

  # Set analysis variables
  chrom_names <- Battenberg::get.chrom.names(imputeinfofile, TRUE)
  logr_file <- paste(tumourname, "_mutantLogR.tab", sep = "")
  externalhaplotypeprefix <- ifelse(!is.na(externalhaplotypefile), paste0(normalname, "_external_haplotypes_chr"), NA)

  # 1. Run haplotying
  #doParallel::registerDoParallel(cores = nthreads)
  #foreach::foreach(i = 1:length(chrom_names), .errorhandling = "remove") %do% {
  for(i in 1:length(chrom_names)) {
    chrom <- chrom_names[i]
    logdebug(paste0("Haplotyping chromosome ", chrom))

    suppressMessages(
      Battenberg::run_haplotyping(
        chrom = chrom,
        tumourname = tumourname[sampleidx],
        normalname = normalname,
        ismale = ismale,
        imputeinfofile = imputeinfofile,
        problemloci = problemloci,
        impute_exe = impute_exe,
        min_normal_depth = min_normal_depth,
        chrom_names = chrom_names,
        snp6_reference_info_file = snp6_reference_info_file,
        heterozygousFilter = heterozygousFilter,
        usebeagle = usebeagle,
        beaglejar = beaglejar,
        beagleref = gsub("CHROMNAME", chrom, beagleref.template),
        beagleplink = gsub("CHROMNAME", chrom, beagleplink.template),
        beaglemaxmem = beaglemaxmem,
        beaglenthreads = beaglenthreads,
        beaglewindow = beaglewindow,
        beagleoverlap = beagleoverlap,
        externalhaplotypeprefix = externalhaplotypeprefix,
        use_previous_imputation = (sampleidx > 1)
      )
    )

    loginfo(paste0("Battenberg RUN HAPLO COMPLETE for ", chrom))
  }
  # doParallel::stopImplicitCluster()

  # 2. Combine all the BAF output into a single file
  Battenberg::combine.baf.files(
    inputfile.prefix = paste(tumourname[sampleidx], "_chr", sep = ""),
    inputfile.postfix = "_heterozygousMutBAFs_haplotyped.txt",
    outputfile = paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
    chr_names = chrom_names
  )

  # Raise error if haplotyping fails to yield tumor BAFs
  tryCatch(
    nrow(fread(paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep = ""))),
    error = function(e) {
      print("Error: Battenberg haplotyping did not yield Tumor BAF results. Quitting.")
      stop(e)
    }
  )

  # 3. Segment the phased and haplotyped BAF data
  Battenberg::segment.baf.phased(
    samplename = tumourname[sampleidx],
    inputfile = paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
    outputfile = paste(tumourname[sampleidx], ".BAFsegmented.txt", sep = ""),
    prior_breakpoints_file = prior_breakpoints_file,
    gamma = segmentation_gamma,
    phasegamma = phasing_gamma,
    kmin = segmentation_kmin,
    phasekmin = phasing_kmin,
    calc_seg_baf_option = calc_seg_baf_option
  )

  # 4. Fit a clonal copy number profile
  Battenberg::fit.copy.number(
    samplename = tumourname[sampleidx],
    outputfile.prefix = paste(tumourname[sampleidx], "_", sep = ""),
    inputfile.baf.segmented = paste(tumourname[sampleidx], ".BAFsegmented.txt", sep = ""),
    inputfile.baf = paste(tumourname[sampleidx], "_mutantBAF.tab", sep = ""),
    inputfile.logr = logr_file,
    dist_choice = clonality_dist_metric,
    ascat_dist_choice = ascat_dist_metric,
    min.ploidy = min_ploidy,
    max.ploidy = max_ploidy,
    min.rho = min_rho,
    min.goodness = min_goodness,
    uninformative_BAF_threshold = uninformative_BAF_threshold,
    gamma_param = platform_gamma,
    use_preset_rho_psi = use_preset_rho_psi,
    preset_rho = preset_rho,
    preset_psi = preset_psi,
    read_depth = 30,
    analysis = "paired"
  )

  # 5. Go over all segments, determine which segements are a mixture of two states and fit a second CN state
  Battenberg::callSubclones(
    sample.name = tumourname[sampleidx],
    baf.segmented.file = paste(tumourname[sampleidx], ".BAFsegmented.txt", sep = ""),
    logr.file = logr_file,
    rho.psi.file = paste(tumourname[sampleidx], "_rho_and_psi.txt", sep = ""),
    output.file = paste(tumourname[sampleidx], "_copynumber.txt", sep = ""),
    output.figures.prefix = paste(tumourname[sampleidx], "_subclones_chr", sep = ""),
    output.gw.figures.prefix = paste(tumourname[sampleidx], "_BattenbergProfile", sep = ""),
    masking_output_file = paste(tumourname[sampleidx], "_segment_masking_details.txt", sep = ""),
    prior_breakpoints_file = prior_breakpoints_file,
    chr_names = chrom_names,
    gamma = platform_gamma,
    segmentation.gamma = NA,
    siglevel = 0.05,
    maxdist = 0.01,
    noperms = 1000,
    calc_seg_baf_option = calc_seg_baf_option
  )

  # 6. If patient is male, get copy number status of ChrX based only on logR segmentation (due to hemizygosity of SNPs)
  if (ismale) {
    if (check_callChrXSubclones(tumourname[sampleidx])) {
      Battenberg::callChrXsubclones(
        tumourname = tumourname[sampleidx],
        X_gamma = 1000,
        X_kmin = 100,
        genomebuild = GENOMEBUILD,
        AR = TRUE,
        chrom_names = chrom_names
      )
    }
  }
}
