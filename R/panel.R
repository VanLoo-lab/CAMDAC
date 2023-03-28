#' Make CAMDAC methylation panel from allele counts
#'  Methylation fractions are obtained by summing M and UM reads across samples
#' @param ac_files Allele count files from CAMDAC
#' @param ac_props Proportions of each sample to use in panel. If NULL, samples are weighted by their
#'  total number of reads, which equals the sum of M and UM counts. If samples are NA, then
#'  proportions are redistributed.
#' @param min_coverage Minimum coverage for a sample's site to be included in panel
#' @param min_samples Minimum number of samples with coverage for a site to be included in panel
#' @param max_sd Maximum standard deviation of methylation for a site to be included in panel
#' @param drop_snps Boolean. If TRUE, drop per-sample CG-SNPs (BAF < 0.1 or BAF > 0.9) from panel
#' @param cores Number of cores to use for calculating HDI
#' @export
panel_meth_from_counts <- function(ac_files, ac_props = NULL, min_coverage = 3, min_samples = 1,
                                   max_sd = 0.1, drop_snps = FALSE, cores = 5) {
  # Load AC files as list, ordering each sample by the same CpG positions
  # Adds a PASS field for us to track and set sites to NA based on filters
  acl <- load_panel_ac_files(ac_files)

  # Apply per-sample CpG constraints
  acl <- apply_coverage_filter(acl, min_coverage)
  acl <- apply_snp_filter(acl, drop_snps)

  # Apply panel CpG constraints
  mask_1 <- min_sample_cg_threshold(acl, min_samples)
  mask_2 <- max_sd_threshold(acl, max_sd)
  panel_mask <- mask_1 & mask_2
  if (sum(panel_mask) == 0) {
    stop("No sites meet panel constraints. Try lowering min_samples or increasing max_sd.")
  }

  # Filter CG sites by panel mask
  acl <- lapply(acl, function(e) e[panel_mask, ])

  # Combine counts to create methylation panel
  panel <- panel_meth_counts(acl, ac_props)

  # Add methylation HDI to panel
  hdi <- calculate_counts_hdi(panel$M, panel$UM, n_cores = cores)
  panel <- cbind(panel, hdi)

  # Return panel object
  return(panel)
}

#' Load allele count files
#' @param ac_files Allele count files from CAMDAC
#' @return List of data tables for each allele counts file
load_panel_ac_files <- function(ac_files) {
  # Set fields to draw from AC files
  ac_load_fields <- c(
    "chrom", "start", "end", "POS", "ref", "alt",
    "total_depth", "M", "UM", "m", "total_counts_m", "BAF"
  )
  # Load ac files as list of data tables
  data <- lapply(ac_files, function(x) {
    v <- data.table::fread(x, select = ac_load_fields)
    setkey(v, chrom, start, end)
    return(v)
  })
  # Find unique cpg positions in all samples so we can create
  #  a shared mapping for the dataset
  uac <- unique_cpg_pos(data)

  # Set all data tables to have the same cpg positions
  cix <- lapply(data, get_overlap_ix, x = uac)
  dix <- lapply(
    seq_along(data),
    function(i) data[[i]][cix[[i]], ]
  )

  # Add PASS field
  dix <- lapply(dix, function(e) {
    e$PASS <- TRUE
    e[is.na(total_counts_m) | is.na(m), PASS := FALSE]
    return(e)
  })

  # Add filename (no path prefixes) to list
  names(dix) <- fs::path_file(ac_files)

  return(dix)
}

unique_cpg_pos <- function(cg_list) {
  x <- lapply(cg_list, function(e) e[, c("chrom", "start", "end")])
  x <- Reduce(rbind, x)
  x <- unique(x)
  return(x)
}

get_overlap_ix <- function(x, y) {
  # Returns row indexes of y that exactly overlap x
  foverlaps(x, y, by.x = c("chrom", "start", "end"), which = T, mult = "first", type = "equal")
}

apply_coverage_filter <- function(e, min_coverage) {
  lapply(
    e,
    function(o) o[PASS == T & total_counts_m < min_coverage, PASS := FALSE]
  )
}

apply_snp_filter <- function(acl, drop_snps) {
  set_snps_na <- function(x) {
    x[
      PASS == T &
        (!is.na(POS) & !is.na(BAF) & dplyr::between(BAF, 0.1, 0.9)),
      PASS := FALSE
    ]
  }
  # If we are to drop SNPS, set all PASS to False where SNPs are called
  if (drop_snps) {
    res <- lapply(acl, set_snps_na)
  } else {
    res <- acl
  }
  return(res)
}

min_sample_cg_threshold <- function(x, min_samples) {
  rs <- Reduce(
    cbind,
    lapply(x, function(o) o$PASS)
  ) %>% rowSums()
  return(rs >= min_samples)
}

max_sd_threshold <- function(x, max_sd) {
  rs <- Reduce(
    cbind,
    lapply(x, function(o) o$m)
  ) %>% rowSds(na.rm = T)
  bool <- ifelse(!is.na(rs) & rs <= max_sd, TRUE, FALSE)
  return(bool)
}

panel_meth_counts <- function(x, ac_props = NULL) {
  # x is a list of data tables with the same cpg positions and fields from CAMDAC ac files
  # mask is a boolean for CpG sites to be included in analysis

  # Set non-passing counts to NA
  x <- lapply(x, function(o) {
    o[PASS == FALSE, `:=`(
      M = NA, UM = NA, m = NA, total_counts_m = NA
    )]
    return(o)
  })

  # If we have been given count proportions, use them to weight methylation rates
  if (!is.null(ac_props)) {
    # Get methylation rates as matrix
    m <- Reduce(
      cbind,
      lapply(seq_along(x), function(i) x[[i]]$m)
    )
    m <- as.matrix(m)

    # Recalculate, weighting by proportions of present data
    pmat = matrix(rep(ac_props, nrow(m)), byrow=T, ncol=length(ac_props))

    # Adjust proportions where beta is NA
    pmat[is.na(m)] = 0
    pmat = pmat/rowSums(pmat)
    m[is.na(m)] = 0 # Allows us to multiply safely with pmat

  # Get new beta based on linear combination of new proportions
  m <- as.numeric(rowSums(m*pmat))

  total_counts_m <- Reduce(
    cbind,
    lapply(seq_along(x), function(i) x[[i]]$total_counts_m) # Get complete counts
  ) %>% rowSums(na.rm = T)

    M <- round(m * total_counts_m, 0)
    UM <- total_counts_m - M
    POS <- NA
    total_depth <- NA
    BAF <- NA
    
    # We can get chrom start and end from one sample as all CpGs aligned before passing to this function
    chrom <- x[[1]]$chrom
    start <- x[[1]]$start
    end <- x[[1]]$end
      } else {
    # Otherwise, sum the counts
    M <- Reduce(cbind, lapply(x, function(o) o$M)) %>% rowSums(na.rm = T)
    UM <- Reduce(cbind, lapply(x, function(o) o$UM)) %>% rowSums(na.rm = T)
    m <- M / (M + UM)
    total_counts_m <- M + UM
    POS <- NA
    total_depth <- NA
    BAF <- NA
    chrom <- x[[1]]$chrom
    start <- x[[1]]$start
    end <- x[[1]]$end
  }

  # Return panel data
  res <- data.table(
    chrom = chrom,
    start = start,
    end = end,
    M = M,
    UM = UM,
    m = m,
    cov = total_counts_m
  )

  return(res)
}

#' Make CAMDAC methylation panel from a matrix of beta values
#' @param mat Matrix of beta values. Rows are CpGs, columns are samples
#' @param chrom Vector of chromosome names
#' @param start Vector of CpG start positions
#' @param end Vector of CpG end positions
#' @param cov Vector of coverage values to give each CpG site. If a matrix is provided, coverage is calculated as the sum of reads for each site.
#' @param cores Number of cores to use for calculating HDI
#' @param min_samples Minimum number of samples that must have a non-NA value for a CpG site to be included in panel
#' @param max_sd Maximum standard deviation of methylation for a site to be included in panel.
#' @export
panel_meth_from_beta <- function(
  mat, chrom , start, end, cov, props, cores, min_samples=1, max_sd=1
){
  # Format chromosome as expected
  chrom = gsub("chr", "", chrom)

  # Get expected formats
  stopifnot(length(props)==ncol(mat))
  mat = as.matrix(mat)

  # Apply min sample filter to CpG sites
  mask_min_samples <- rowSums(!is.na(mat))>=min_samples
  mat = mat[mask_min_samples,]

  # Apply max sd filter to CpG sites
  mask_max_sd <- rowSds(mat, na.rm=T) <= max_sd
  mask_max_sd[is.na(mask_max_sd)] <- TRUE
  mat = mat[mask_max_sd,]

  # Apply filter to coverage depending on whether it is a single value, vector or matrix
  if(is.null(dim(cov))){

    if(length(cov)==1){
      pass
    }else{
      cov = cov[ mask_min_samples | mask_max_sd]
    }

  } else {
    cov = cov[ mask_min_samples | mask_max_sd,]
    cov = rowMeans(cov, na.rm=T)
  }

  # Set proportions as matrix
  pmat = matrix(rep(props, nrow(mat)), byrow=T, ncol=length(props))

  # Adjust proportions where beta is NA
  pmat[is.na(mat)] = 0
  pmat = pmat/rowSums(pmat)
  mat[is.na(mat)] = 0 # Allows us to multiply safely with pmat

  # Get new beta based on linear combination of new proportions
  nbeta <- as.numeric(rowSums(mat*pmat))
  # Get new counts based on coverage
  M <- round(cov*nbeta, 0)
  UM <- cov-M
  
  # Set panel
  panel <- data.table(
    chrom = chrom,
    start = start,
    end = end,
    M = M,
    UM = UM,
    m = nbeta,
    cov = cov
  )
  # Add methylation HDI to panel
  hdi <- calculate_counts_hdi(panel$M, panel$UM, n_cores = cores)
  panel <- cbind(panel, hdi)

  return(panel)
}
