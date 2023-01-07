#' Pre-process CAMDAC methylation data

# Working knowledge
## In the below exampls, PX is the patient ID and R1 is the tumour region
# dt_sample_v_normal_m.RData - No longer created in latest version, in favour of dt_tumour_and_normal_m.
# dt_tumour_m.RData - Base bulk methylation data: with columns: chromosome,start,end,M_b,UM_b,m_b,cov_b,SNP_b.
#   No longer output in latest CAMDAC version (run_methylation_data...*.R)
# dt_normal_m.RData - This has the same data as the dt_tumour_m but includes HDI for later merging.
#   I keep things consistent by including.
# purified_tumour_DMP_calls.RData
# purified_tumour.RData
# PX1.R1.DMPs.txt - summary stats for the number of DMPs and tumour-normal differential methylation data
# PX1_R1_methylation_rate_summary.pdf
# Rplots.pdf

#

# Code ----

#' Calculate HDI interval width
#' @noRd
intervalWidth <- function(lowTailPr, ICDFname, credMass, ...) {
  ICDFname(credMass + lowTailPr, ...) - ICDFname(lowTailPr, ...)
}

#' HDI of ICDF
#' @param ICDFname The inverse cumulative density function of the distribution.
#' @param credMass The desired mass of the HDI region.
#' @param tol Tolerance parameter for optimisation. the lower the tolerance,the
#'   longer the optimisation, but the higher the accuracy.
#'   According to CAMDAC RRBS comments, tol=1e-4 gives values
#'   of the same accuracy as our max resolution.
#'   This function is adapted from Greg Snow's TeachingDemos package
#'   E.g.Determine HDI of a M=30 and UM=12 CpG
#'   Adding 1 to shape parameter ensures uniform beta(1,1) is updated with our counts
#'   HDIofICDF(qbeta,shape1 = 30+1 , shape2 = 12+1 )
#' @return Highest density interval (HDI) limits in a vector.
#' @keywords internal
HDIofICDF <- function(ICDFname, credMass = 0.99, tol = 1e-4, ...) {
  incredMass <- 1.0 - credMass

  # Here, shape parameters are passed to ICDFname function via `...`
  optInfo <- optimize(f = intervalWidth, interval = c(0, incredMass), ICDFname = ICDFname, credMass = credMass, tol = tol, ...)

  HDIlowTailPr <- optInfo$minimum
  vec <- setNames(object = ICDFname(c(HDIlowTailPr, credMass + HDIlowTailPr), ...), nm = c("low", "high"))
  return(data.frame(lo = vec[[1]], high = vec[[2]]))
  # return(vec)
}

# Calculate HDI counts for unique combinations of records to speed up processing time
unique_calculate_counts_hdi <- function(M, UM, n_cores = 1, itersplit = 5e5) {
  # Itersplit default: Benchmarking found 500K cpgs takes ~1min
  # Validate M-length
  inp_len <- length(M)
  stopifnot(inp_len == length(UM))

  # Get unique pairs to save computation
  udata <- unique(data.table(M = M, UM = UM))
  M <- udata$M
  UM <- udata$UM
  inp_len <- length(M)
  rm(udata)

  # Split data for parallel chunks
  split_factor <- make_split_factor(inp_len, itersplit)
  M <- split(M, split_factor)
  UM <- split(UM, split_factor)

  # Calculate HDI parallel
  doParallel::registerDoParallel(cores = n_cores)
  # mapply is used to vectorise over M and UM, which are arrays
  hdi <- foreach(M = M, UM = UM, .combine = "c") %dopar% {
    mapply(
      FUN = HDIofICDF,
      shape1 = M + 1,
      shape2 = UM + 1,
      MoreArgs = list(ICDFname = qbeta),
      SIMPLIFY = FALSE,
      USE.NAMES = TRUE
    )
  }
  doParallel::stopImplicitCluster()

  # Bind result as data.table
  hdi <- data.table::rbindlist(hdi)
  res <- cbind(
    data.table(M = unlist(M), UM = unlist(UM)),
    hdi
  )

  return(res)
}

calculate_counts_hdi <- function(M, UM, n_cores = 1, itersplit = 5e5) {
  # Calculate HDI and bind to original data. Adds columns "m_x_low" and "m_x_high"
  loginfo("Calculating HDI from counts.")
  u_hdi <- unique_calculate_counts_hdi(M, UM, n_cores = n_cores, itersplit = itersplit)
  u_hdi <- round(u_hdi, digits = 5)
  # Combine original data with HDI in order
  loginfo("Merging HDI result")
  hdi_data <- merge(
    data.table(M = M, UM = UM, i = seq(length(M))),
    u_hdi,
    all.x = TRUE,
    by = c("M", "UM")
  )
  names(hdi_data) <- c("M", "UM", "hdi_i", "m_x_low", "m_x_high")
  hdi_data <- hdi_data[order(hdi_i), .(m_x_low, m_x_high)]
  # Return HDI data alone (able to cbind original data)
  return(hdi_data)
}

# Return allele counts data restricted to methylation sites and formatted for
# downstream deconvolution and DMP identification
process_methylation <- function(allele_counts, min_meth_loci_reads = 3) {
  # Limit data to CpG/CCGG sites with methylation data.
  methyl <- allele_counts[width > 1 & !is.na(m) & total_counts_m > min_meth_loci_reads]
  rm(allele_counts)

  # Annotate heterozygous SNPs for downstream CG-SNP investigation
  # FEATURE: BAF is already counted in allele_counts,
  # therefore this label may be more appropriate earlier in the pipeline
  methyl[, SNP := fifelse(!is.na(BAF) & (BAF >= 0.15) & (BAF <= 0.85), 1, 0)]

  # Add column of methylation coverage as cov. #TODO: More efficient DT rename here?
  methyl[, cov := total_counts_m]

  # Add CG-SNP status (CG-forming or destroying). Required for accurate CG-copy number assignment
  methyl[, cg_snp := classify_cg_snp(start, width, POS, ref, alt)]

  # Subset to required output columns
  methyl <- methyl[, .(chrom, start, end, M, UM, m, cov, SNP, BAF, cg_snp)]

  # Return data
  return(methyl)
}



save_methylation_df <- function(methyl, sample, config) {
  output_file <- build_output_name(sample, config, "methylation")
  data.table::fwrite(methyl, file = output_file)
}

# A key result of run_methylation_data_processing is the dt_tumour_and_normal_m.RData object
# I would like to see how far I can get without combining these two, but rather working with the data separately (memory issues afterall)
# I.e. can differential methylation analysis start by loading each separately?
combine_tumour_normal_methylation <- function(t_meth, n_meth) {
  new_names <- sapply(names(n_meth), function(x) {
    ifelse(x %in% c("CHR", "chrom", "start", "end"), x, paste0(x, "_n"))
  })
  names(n_meth) <- new_names

  # Set keys for merge (sorts table internally)
  #   We first correct orderings for chrom fields by making them factors
  t_meth$chrom <- factor(t_meth$chrom, levels = c(1:22, "X", "Y"))
  n_meth$chrom <- factor(n_meth$chrom, levels = c(1:22, "X", "Y"))
  setkey(t_meth, chrom, start, end)
  setkey(n_meth, chrom, start, end)

  # Combine into a single CpG table.
  #   All normal fields are now prefixed with 'i.'
  #   nomatch=0 drops mismatching fields, rather than retaining them as NA
  #   type="equal" searches for exact CpG range matches
  meth_c <- foverlaps(n_meth, t_meth, nomatch = 0, type = "equal")
  meth_c$i.start <- NULL
  meth_c$i.end <- NULL
  setkey(meth_c, chrom, start, end)
  return(meth_c)
}

annotate_cgs_with_cnas <- function(meth_c, tumour) {
  cna <- tumour$cna
  cna$ascna$chrom <- factor(cna$ascna$chrom, levels = c(1:22, "X", "Y"))
  setkey(cna$ascna, chrom, start, end)
  # Add additional columns so segments can be referenced elswhere in code
  cna$ascna$seg_start <- cna$ascna$start
  cna$ascna$seg_end <- cna$ascna$end
  meth_cna <- foverlaps(cna$ascna, meth_c, nomatch = 0)
  meth_cna[, i.start := NULL]
  meth_cna[, i.end := NULL]

  # Set CG copy number us BAF threshold of 0.5.
  #  At CG-SNP sites, reads will only contain CGs depending on whether the site is a CG-forming or CG-destroying SNP.
  #  Therefore the methylation copy number at these loci will reflect the major or minor allele.
  meth_cna[, CG_CN := data.table::fcase(
    #  If not a CG-SNP, take total copy number
    is.na(BAF), CN,
    #  Take major if majority allele contains CG (CG-destroying with low BAF or CG-forming with high BAF)
    (BAF <= 0.5 & cg_snp == "D") | (BAF >= 0.5 & cg_snp == "F"), nA,
    #  Take minor if majority allele does no contain CG (CG-destroying with high BAF or CG-forming with low BAF)
    (BAF > 0.5 & cg_snp == "D") | (BAF < 0.5 & cg_snp == "F"), nB,
    # Set to total CN for any other case
    rep(TRUE, nrow(meth_cna)), CN
  )]

  # Set normal CG copy number.
  # CG_CN_n is 1 at CG-destroying/forming hetrozygous SNPs
  meth_cna[, CG_CN_n := data.table::fifelse(is.na(BAF), 2, 1, na = NA)]
  # Set normal copy number on sex chromosome X in MALES
  # In males, it has CN=1 outside PAR regions and 2 within.
  # TODO: meth_cna <- calculate_cg_cn_norm(meth_cna, patient_sex, reference_genome)

  # Add overall tumour purity
  meth_cna[, p := cna$purity]
  return(meth_cna)
}

# TODO: Include CCGG sites in this function for RRBS
classify_cg_snp <- function(start, width, POS, ref, alt) {
  cg_snp_class <- data.table::fcase(
    # A CG-SNP is any position with a width >1 and a non-na POS
    # Return NA for sites that do not meet this criteria
    width != 2 | is.na(POS), NA_character_,
    # CG-destroying SNPs have a reference C at CG-start or G at CG-end.
    (start == POS & ref == "C") | (start + 1 == POS & ref == "G"), "D",
    # CG-forming SNPs have an alt C at the CG-start or G at the CG-end.
    (start == POS & alt == "C") | (start + 1 == POS & alt == "G"), "F",
    default = NA_character_
  )
  return(
    factor(cg_snp_class, levels = c("F", "D", NA_character_))
  )
}

calculate_mt <- function(mb, mn, p, CN) {
  tumour_frac <- p * CN
  normal_frac <- (1 - p) * 2 # Normal CN assumed to be 2
  mt <- ((
    mb * (tumour_frac + normal_frac)
  ) - (
    mn * normal_frac
  )) /
    tumour_frac
  return(mt)
}

calculate_mt_cov <- function(cov_b, p, CN) {
  # Effective tumour coverage estimated by deconvolving bulk coverage
  # The fractin of reads reporting the tumour is a function of the tumour purity,
  #   and copy number. If CN was not included two sites with the same purity, cov_b and CN_norm can
  #   differ by CN and be incorrectly inferred
  tumour_frac <- p * CN
  normal_frac <- (1 - p) * 2 # Normal CN assumed to be 2
  cov_t <- round(
    cov_b * tumour_frac / (tumour_frac + normal_frac),
    digits = 0
  )
  return(cov_t)
}

deconvolve_bulk_methylation <- function(meth_c) {
  # Deconvolve methylation
  meth_c[, m_t_raw := calculate_mt(
    m, m_n, p, CG_CN
  )]

  # Correct pure tumour methylation rates set outside 0 and 1 after deconvolution
  meth_c[, m_t := data.table::fcase(
    m_t_raw < 0, 0,
    m_t_raw > 1, 1,
    rep(TRUE, nrow(meth_c)), m_t_raw
  )]

  # Calculate tumour coverage by deconvolution
  meth_c[, cov_t := calculate_mt_cov(cov, p, CG_CN)]

  return(meth_c)
}

filter_deconvolved_methylation <- function(meth_c) {
  meth_c[
    CN > 0 & # Remove homozygous deletions
      cov_t >= 3 & # Remove low-coverage CpGs (after deconvolution)
      !is.na(m_t_raw) # Capture any errors in deconvolution. Should be none!
  ]
}

#' Calculate HDI by simulation
#' @keywords internal
calculate_m_t_hdi <- function(meth_c, n_cores, itersplit = 1e5) {
  inp_len <- nrow(meth_c)
  split_factor <- make_split_factor(inp_len, itersplit)

  msplit <- iterators::isplit(meth_c, split_factor)

  # Calculate HDI
  doParallel::registerDoParallel(cores = n_cores)
  hdi_all <- foreach(v = msplit, .combine = "rbind") %dopar% {
    x <- v$value
    hdi <- vec_HDIofMCMC_mt(
      M_b = x$M,
      UM_b = x$UM,
      M_n = x$M_n,
      UM_n = x$UM_n,
      p = x$p,
      CN = x$CG_CN,
      credMass = 0.99
    )
    colnames(hdi) <- c("m_t_low", "m_t_high")
    return(hdi)
  }
  doParallel::stopImplicitCluster()

  meth_c <- cbind(meth_c, hdi_all)
  return(meth_c)
}

#' Calculate HDI by simulation
#'
#' Computes highest density interval from a sample of representative values,
#'  estimated as shortest credible interval for a unimodal distribution
#'
#' @param M_b counts methylated in the tumour
#' @param UM_b counts unmethylated in the tumour
#' @param M_n counts methylated in the normal
#' @param UM_n counts unmethylated in the normal
#' @param p tumour purity
#' @param CN total tumour copy number
#' @param CN_n total normal copy number
#' @param credMass default is 0.99
#' credMass is a scalar between 0 and 1, indicating the mass within the
#' credible interval that is to be estimated.
#' @return Value: HDIlim is a vector containing the limits of the HDI
#' @keywords internal
HDIofMCMC_mt <- function(M_b, UM_b, M_n, UM_n, p, CN, credMass = 0.99) {
  # Simulate beta distributions from bulk and normal methylation
  bulk_dist <- rbeta(n = 2000, shape1 = M_b + 1, shape2 = UM_b + 1)
  normal_dist <- rbeta(n = 2000, shape1 = M_n + 1, shape2 = UM_n + 1)

  # Get constants for effect of purity and copy number
  tumour_frac <- p * CN
  normal_frac <- (1 - p) * 2 # Normal CN assumed to be 2
  bulk_constant <- tumour_frac + normal_frac

  # Deconvolve methylation rates from monte carlo simulation to approximate mt distribution
  mt_dist <- ((bulk_dist * bulk_constant) - (normal_dist * normal_frac)) / tumour_frac

  # Calculate the HDI of mt from the simulation
  # See (Kruschke, K., 2015, Doing Bayesian Data Analysis, 721–736)
  mt_dist <- sort(mt_dist)
  # get the width of the nth percentile where n=credMass
  ci_idx_length <- ceiling(credMass * length(mt_dist))
  # get the diff between all pairs at a suitable width
  ciWidth <- diff(x = mt_dist, lag = ci_idx_length)
  # Return HDI
  HDIlim <- (c(
    lower = mt_dist[which.min(ciWidth)],
    upper = mt_dist[which.min(ciWidth) + ci_idx_length]
  ))

  return(HDIlim)
}

# Vectorise HDI simulation
vec_HDIofMCMC_mt <- function(...) {
  vf <- Vectorize(HDIofMCMC_mt) # Accepts vectors with dim(y,0)
  result <- vf(...) # Returns a matrix with dim(x ,y)
  # Transpose result to return rows corresponding to original data: dim(y,x)
  return(t(result))
}

# TODO: Calculate HDI and annotate data table.
# Note: CAMDAC defines a column name suffix depending on the tumour/normal status

# MISC DEV -----


# SETUP TEST ( CAMP ) Nov 2021:
# devtools::load_all("~/projects/nmns93-camdac")
# opt = load_camdac_opts_from_input("DTB-042", "~/data/211027_CAMDAC/211001_camdac_wgbs_all_inputs.tsv",
#                                   "~/data/211027_CAMDAC/", "~/projects/nmns93-camdac/tests/testthat/data/pipeline_files/")
# tumour=opt[[1]]; normal=opt[[2]]; config=opt[[3]]
# ac_file = build_output_name(tumour, config, "allele_counts")
# allele_counts = data.table::fread(ac_file, nrow=1e5)
# min_meth_loci_reads=3
# t_m = data.table::fread(build_output_name(tumour, config, "methylation"))
# t_n=data.table::fread(build_output_name(normal, config, "methylation"))
#  qplot(x=m, geom="density", data=t_n) + geom_density(data=t_m, mapping=aes(x=m), colour="pink")

# SETUP TEST (LOCAL PC):
# devtools::load_all("~/Documents/CAMDAC")
# setwd("~/Documents/CAMDAC/tests/testthat")
# Note files will soon change to csv.gz
# allele_counts = readRDS("runs/pgp-test/Allelecounts/pgp-001/pgp-test.pgp-001.SNPs.CpGs.all.sorted.rds")
# allele_counts = add_fake(allele_counts)
# min_meth_loci_reads=3

# TEMP: ADD fake data for testing while full dataframe is insufficient
add_fake <- function(data) {
  # TODO: As there are no SNP-cpgs in my test data (0.25% occurence in full data from RRBS inspection)
  #    ensure that add_fake injection is migrated to test file.
  fake_sites <- list(
    # Add sites for cpg and snps co-occuring
    list("chr1", 1, 860720, 860721, 2, 860720, "C", "T", 0, 7, 7, 0, 11, 0, 11, 5, 6, 11, 5 / 11, 0, 0, 0, 7, 4, 0, 0, 0, 0, 0, 0, 0, 0, "30,44,60,60,30,30,60,6,60,60,60"),
    list("chr10", 10, 126707095, 126707096, 2, 126707095, "C", "T", 0, 7, 7, 0, 11, 0, 11, 5, 6, 11, 5 / 11, 0, 0, 0, 7, 4, 0, 0, 0, 0, 0, 0, 0, 0, "30,44,60,60,30,30,60,6,60,60,60"),
    list("chr2", 2, 241489488, 241489489, 2, 241489488, "C", "T", 0, 7, 7, 0, 11, 0, 11, 5, 6, 11, 5 / 11, 0, 0, 0, 7, 4, 0, 0, 0, 0, 0, 0, 0, 0, "30,44,60,60,30,30,60,6,60,60,60"),

    # Add sites to test loci where hits occur on both cpgs
    list("chr7", 7, 1652650, 1652651, 2, 1652650, "C", "T", 0, 7, 7, 0, 11, 0, 11, 5, 6, 11, 5 / 11, 0, 0, 0, 7, 4, 0, 0, 0, 0, 0, 0, 0, 0, "30,44,60,60,30,30,60,6,60,60,60"),
    list("chr7", 7, 1652650, 1652651, 2, 1652651, "C", "T", 0, 7, 7, 0, 11, 0, 11, 5, 6, 11, 5 / 11, 0, 0, 0, 7, 4, 0, 0, 0, 0, 0, 0, 0, 0, "30,44,60,60,30,30,60,6,60,60,60")
  )
  sites <- rbindlist(fake_sites)
  names(sites) <- names(data)
  rbind(data, sites)
}

# TODO: Upstream, perhaps in filter_bad_ac_rows (maybe turn into ac_quality_filter).
# Remove super-high coverage sites. Lili uses 0.9999 quantile,
# and agreed in lab meeting this is acceptabel.

# TODO: Remove sites duplicated due to two co-localising SNP occurring at both CG nucleotide loci
# OR CG and CCGG.
# Note: CRRBS deduplicates here and recalculates the methylation rates. I would instead
# ensure this SNP protocol is handled at allele-counts section.


# This was raising roxygen error. We have this in allele_counts so not sure why it's here
# #' @param dt Data table with chrom, start and end columns
# drop_all_duplicates <- function(dt){
#
#   # Remove sites duplicated due to multiple SNPs overlapping CG/CCGG or
#   # SNPs on both CG nucleotides.
#   # FUTURE: CAMDAC to allow a one-to-many relationship between CGs and SNPs
#   # pileup_summary <- drop_all_duplicates(pileup_summary)
#
#   # We want to remove all CG records that occur more than once, but the `duplicate` function doesn't return the first
#   # of the duplicated records. Here, we use the fromLast argument to search for duplicates in the opposite direction,
#   # and merge these results with `|` to capture all duplicate instances.
#
#   dt [ !(
#     duplicated(dt, by=c("chrom","start","end"), fromLast=F) |
#       duplicated(dt, by=c("chrom","start","end"), fromLast=T)
#   )]
#
# }

# Helper function to split data
make_split_factor <- function(nrows, itersplit) {
  # Get the integer by which we will split data
  # If nrows is < itersplit, set the factor to 1
  split_factor <- ifelse(
    nrows < itersplit,
    1,
    round(nrows / itersplit, 0)
  )

  # Repeat a sequence of our split factor
  # then sort it to ensure data is split in order
  split_f <- sort(
    rep_len(
      seq(split_factor),
      nrows
    )
  )

  return(split_f)
}
