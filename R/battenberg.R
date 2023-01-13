
# This is parallelised per chromomse in bb, so will run over multiple files in the future
allele_counts_to_beagle_vcf_ <- function(allele_counts, outfile, impute_legend, min_het_BAF = 0.15) {
  stopifnot(min_het_BAF < 0.5) # Ensure threshold is the lower of the BAF cutoff

  # Read allele counts file
  normal_counts <- fread_chrom(allele_counts)

  # Filter for impute legend SNPs
  known_snps <- data.table::fread(impute_legend)[, 1:4]
  stopifnot(names(known_snps) == c("id", "position", "a0", "a1"))
  data.table::setnames(known_snps, "position", "POS")
  normal_counts <- normal_counts[known_snps, , on = "POS"][!is.na(Good_depth)]
  data.table::setnames(normal_counts, "a0", "ref")
  data.table::setnames(normal_counts, "a1", "alt")
  nucleotides <- c("A", "T", "C", "G")
  normal_counts <- normal_counts[ref %in% nucleotides & alt %in% nucleotides] # Remove multi-allelic
  normal_counts <- normal_counts[Good_depth > 0]

  # Calculate BAF
  normal_counts[, ref_names := paste0("Count_", ref)]
  normal_counts[, alt_names := paste0("Count_", alt)]
  normal_counts[, ref_counts := data.table::fcase(ref_names == "Count_A", Count_A, ref_names == "Count_C", Count_C, ref_names == "Count_G", Count_G, ref_names == "Count_T", Count_T)]
  normal_counts[, alt_counts := data.table::fcase(alt_names == "Count_A", Count_A, alt_names == "Count_C", Count_C, alt_names == "Count_G", Count_G, alt_names == "Count_T", Count_T)]
  normal_counts[, BAF := alt_counts / (alt_counts + ref_counts)]
  normal_counts <- normal_counts[!is.na(BAF)]

  # Assign VCF genotype string
  normal_counts[, SAMP001 := data.table::fcase(
    BAF <= min_het_BAF, "0/0",
    BAF >= 1 - min_het_BAF, "1/1",
    BAF > min_het_BAF & BAF < 1 - min_het_BAF, "0/1"
  )]

  # Format VCF for writing
  data.table::setnames(normal_counts, "#CHR", "#CHROM")
  data.table::setnames(normal_counts, "ref", "REF")
  data.table::setnames(normal_counts, "alt", "ALT")
  normal_counts[, ID := "."]
  normal_counts[, QUAL := "."]
  normal_counts[, FILTER := "PASS"]
  normal_counts[, INFO := "."]
  normal_counts[, FORMAT := "GT"]

  vcf <- normal_counts[, c("#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT", "SAMP001")]

  # vcf = lift_hg19_to_hg38(vcf)
  # Write beagle input VCF
  writevcf.beagle(vcf, outfile)
}

# Convert CAMDAC allele counts object to per-chromosome allele frequency files for battenberg
create_allele_frequencies_chr <- function(camdac_allele_counts, outdir, sample_af_prefix, min_depth) {
  # Load data table of SNP positions
  dt <- fread_chrom(camdac_allele_counts, select = c("chrom", "POS", "ref", "alt", "ref_counts", "alt_counts"))[!is.na(POS)]
  # Filter for min depth
  dt <- dt[ref_counts + alt_counts >= min_depth]
  # Split per chromosome for parallel processing
  dt_per_chrom <- split(dt, dt$chrom)
  rm(dt)

  af_files <- foreach(dt_chrom = dt_per_chrom) %dopar% {
    chromosome <- unique(dt_chrom$chrom)[[1]]
    # Set output file names with correct suffix for X and Y chromosomes
    outfile_chrom <- data.table::fcase(
      chromosome == "X", "23",
      chromosome == "Y", "24",
      !(chromosome %in% c("X", "Y")), chromosome
    )
    outfile <- fs::path(outdir, sprintf("%s_alleleFrequencies_chr%s.txt", sample_af_prefix, outfile_chrom))

    # Set CHR column to allelefreq format expected
    data.table::setnames(dt_chrom, "chrom", "#CHR")

    # Label alleles to allelefreq hedings. E.g. T -> Count_T
    dt_chrom[, `:=`(ref = paste0("Count_", ref), alt = paste0("Count_", alt))]

    # Pivot ref labels as columns with counts per SNP position as values
    # Note: We perform for ref and alt separately, so that we can later stack and sum the count columns by SNP position
    cubed <- data.table::cube(dt_chrom, .(result = sum(ref_counts)), by = c("#CHR", "POS", "ref"))
    ref_table <- data.table::dcast(cubed, `#CHR` + POS ~ ref, value.var = "result")[!is.na(`#CHR`) & !is.na(POS)]
    ref_table[is.na(ref_table)] <- 0

    # Pivot alt labels as columns with counts per SNP position as values
    cubed <- data.table::cube(dt_chrom, .(result = sum(alt_counts)), by = c("#CHR", "POS", "alt"))
    alt_table <- data.table::dcast(cubed, `#CHR` + POS ~ alt, value.var = "result")[!is.na(`#CHR`) & !is.na(POS)]
    alt_table[is.na(alt_table)] <- 0

    # If data is sparse, some SNPs may not be present
    # fix missing fields in ref and alt table
    fix_count_fields <- function(data) {
      count_fields <- paste0("Count_", c("A", "C", "G", "T"))
      missing <- setdiff(count_fields, names(data))
      if (length(missing) > 0) {
        data[, missing] <- 0
      }
      return(data)
    }
    ref_table <- fix_count_fields(ref_table)
    alt_table <- fix_count_fields(alt_table)

    # Stack ref and alt count columns
    headings <- c("#CHR", "POS", "Count_A", "Count_C", "Count_G", "Count_T", "Good_depth")
    headings_sub <- headings[1:6] # Good_depth not implemented. Calculated downstream as sum
    allele_counts <- rbind(ref_table[, ..headings_sub], alt_table[, ..headings_sub])
    rm(dt_chrom, cubed, ref_table, alt_table)

    # Sum Count_Nucleotide columns by SNP position and calculate Good_depth field
    allele_counts <- allele_counts[, .(Count_A = sum(Count_A), Count_C = sum(Count_C), Count_G = sum(Count_G), Count_T = sum(Count_T)), by = .(`#CHR`, POS)]
    allele_counts[, Good_depth := sum(Count_A, Count_C, Count_G, Count_T), by = 1:nrow(allele_counts)]
    allele_counts <- allele_counts[order(POS)] # Order by POS only as only a single chrom is used.

    # Write to output file
    write.table(allele_counts, file = outfile, sep = "\t", row.names = F, col.names = T, quote = F)
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

    tumour_dt <- fread_chrom(tumour_af)
    setkey(tumour_dt, "#CHR", "POS")
    normal_dt <- fread_chrom(normal_af)
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

camdac_to_battenberg_allele_freqs <- function(camdac_tumour_ac, tumour_prefix, camdac_normal_ac, normal_prefix, outdir, min_normal_depth) {
  create_allele_frequencies_chr(camdac_tumour_ac, outdir, tumour_prefix, min_depth = 0) # TODO make function(exists)
  create_allele_frequencies_chr(camdac_normal_ac, outdir, normal_prefix, min_depth = min_normal_depth)
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

  # TODO: Create the alleleCounts.tab file? Not necessary for pipeline.

  # TODO: Create ASCAT plots and simplify the tables by passing sample name to the write output function:
  # see: https://github.com/Wedge-lab/battenberg/blob/43686673566cf5adbd8d00e2450d70eced27696d/R/prepare_wgs.R#L154

  return(filenames)
}

check_callChrXSubclones <- function(TUMOURNAME) {
  # Helper function to check whether we can run callChrX subclones i.e. enough SNPs in non-par regions
  PCFinput <- data.frame(read_table_generic(paste0(TUMOURNAME, "_mutantLogR_gcCorrected.tab")), stringsAsFactors = F)
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
                                    nthreads = 1,
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
  doParallel::registerDoParallel(cores = nthreads)
  foreach::foreach(i = 1:length(chrom_names), .errorhandling = "remove") %dopar% {
    chrom <- chrom_names[i]
    log_debug(paste0("Haplotyping chromosome ", chrom))

    # TODO: Run but suppress std-err from beagle. Use tryCatch to print stderr if error occurs.
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

    print(paste0("Battenberg RUN HAPLO COMPLETE for ", chrom))
  }
  doParallel::stopImplicitCluster()

  # 2. Combine all the BAF output into a single file
  Battenberg::combine.baf.files(
    inputfile.prefix = paste(tumourname[sampleidx], "_chr", sep = ""),
    inputfile.postfix = "_heterozygousMutBAFs_haplotyped.txt",
    outputfile = paste(tumourname[sampleidx], "_heterozygousMutBAFs_haplotyped.txt", sep = ""),
    chr_names = chrom_names
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
    output.file = paste(tumourname[sampleidx], "_subclones.txt", sep = ""),
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
        TUMOURNAME = tumourname[sampleidx],
        X_GAMMA = 1000,
        X_KMIN = 100,
        GENOMEBUILD = GENOMEBUILD,
        AR = TRUE
      )
    }
  }
}
