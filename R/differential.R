# Evan Miller's closed form solution for the probability that
# a draw from a beta dist is greater than another
# Takes counts and methylation fractions for normal and bulk
prob_diff_meth <- function(M_n, UM_n, M, UM) {
  # Return NA when counts are not given
  if (any(is.na(c(M_n, UM_n, M, UM)))) {
    return(NA)
  }

  M_n <- M_n + 1
  UM_n <- UM_n + 1
  M <- M + 1
  UM <- UM + 1
  j <- seq.int(0, round(M) - 1)
  log_vals <- (lbeta(M_n + j, UM_n + UM) - log(UM + j) -
    lbeta(1 + j, UM) - lbeta(M_n, UM_n))
  1 - sum(exp(log_vals))
}
# Vectorized
v_prob_diff_meth <- Vectorize(prob_diff_meth)

# Calculate probability of DMP (difference between betas) from bulk and normal counts
calc_prob_dmp <- function(M_n, UM_n, M, UM, itersplit = 5e5, ncores = 5) {
  split_factor <- make_split_factor(length(M_n), itersplit)
  msplit <- iterators::isplit(seq(length(M_n)), split_factor)

  doParallel::registerDoParallel(cores = ncores)
  prob <- foreach(v = msplit, .combine = "c") %do% {
    x <- v$value
    ph <- v_prob_diff_meth(M_n[x], UM_n[x], M[x], UM[x])
    return(ph)
  }
  doParallel::stopImplicitCluster()

  # Format values
  prob <- data.table::fcase(
    is.na(prob), 0.5,
    prob > 1, 1,
    prob < 0, 0,
    rep_len(TRUE, length(prob)), prob # Otherwise return value
  )

  return(prob)
}

#' Call differentially methylated positions
#' @keywords internal
call_dmps <- function(pmeth, nmeth, effect_size = 0.2, prob = 0.99, itersplit = 5e5, ncores = 5) {
  stopifnot(nrow(pmeth) == nrow(nmeth))

  # Set variables
  M_n <- nmeth$M
  UM_n <- nmeth$UM
  M <- pmeth$M
  UM <- pmeth$UM
  m <- pmeth$m
  m_b_diff <- pmeth$m - nmeth$m
  m_t_diff <- pmeth$m_t - nmeth$m

  phypo <- calc_prob_dmp(M_n, UM_n, M, UM, ncores = ncores, itersplit = itersplit)

  prob_DMP <- data.table::fifelse(m_t_diff > 0, 1 - phypo, phypo) # I.e. if bulk is greater than normal then it's a hyper DMP
  rm(phypo)

  DMP_b <- prob_to_call(prob_DMP, m_b_diff, effect_size = effect_size, prob = prob)
  DMP_t <- prob_to_call(prob_DMP, m_t_diff, effect_size = effect_size, prob = prob)


  res <- cbind(
    pmeth,
    data.table(prob_DMP, DMP_b, DMP_t,
      ndmp_m = nmeth$m,
      ndmp_cov = nmeth$cov,
      ndmp_ml = nmeth$m_x_low,
      ndmp_mh = nmeth$m_x_high
    )
  )

  return(res)
}

prob_to_call <- function(p, mdiff, effect_size = 0.2, prob = 0.99) {
  data.table::fcase(
    p >= prob & mdiff >= effect_size, "hyper",
    p >= prob & mdiff <= (-effect_size), "hypo"
  )
}

#' Add CAMDAC region annotations to dt.
#' DT must have chrom, start, end
#' @noRd
annotate_dmp_regions <- function(dt, all_regions_anno) {
  # Ensure chromosomes are correct format
  dt[, chrom := factor(chrom, levels = c(1:22, "X", "Y"), ordered = TRUE)]
  all_regions_anno[, chrom := factor(chrom, levels = c(1:22, "X", "Y"), ordered = TRUE)]

  # Overlap annotated regions and CpG methylation objects
  setkey(dt, chrom, start, end)
  setkey(all_regions_anno, chrom, start, end)
  dt <- foverlaps(all_regions_anno, dt, type = "any", nomatch = NULL)

  # Order regions
  dt <- dt[order(
    as.numeric(cluster_id),
    factor(chrom, levels = c(1:22, "X", "Y"), ordered = TRUE),
    start, end
  ), ]
  return(dt)
}

#' Count CpGs within DMP annotations
#' @keywords internal
get_cluster_counts <- function(dt) {
  cluster_stats <- dt[, .SD, .SDcols = c("m_t", "m_n", "DMP_t", "cluster_id")]
  cluster_counts <- cluster_stats[, .(
    CpG_counts = length(DMP_t),
    DMP_counts = length(DMP_t[!is.na(DMP_t)]),
    consec_DMPs = max_consec_dmp(DMP_t)
  ), by = "cluster_id"]
  return(cluster_counts)
}

#' Summarise CG stats per DMR
#' @keywords internal
collapse_cpg_to_dmr <- function(dt) {
  dt <- dt[!is.na(DMR) & !is.na(DMP_t)]
  dt <- dt[,
    .(
      m_n = mean(m_n),
      m_n_low = mean(m_x_low_n),
      m_n_high = mean(m_x_high_n),
      m_t = mean(m_t),
      m_t_low = mean(m_t_low),
      m_t_high = mean(m_t_high),
      prob = mean(prob_DMP),
      CG_CN = mean(CG_CN),
      nA = mean(nA, na.rm = T),
      nB = mean(nB, na.rm = T),
      # segment = paste(unique(segment), collapse = ";"),
      DMR_type = set_dmr_type(DMP_t)
    ),
    by = .(cluster_id, chrom, i.start, i.end)
  ]
  setnames(dt, "i.start", "start")
  setnames(dt, "i.end", "end")
  dt[, DMR := "DMR"]
  return(dt)
}

# Reapply DMR annotations, which are currently lost when collapsing CpGs
re_annotate_dmrs <- function(dt, all_regions_anno) {
  all_regions_anno[, `:=`(
    chrom = NULL, start = NULL, end = NULL, CpG_counts = NULL
  )]

  dt <- merge(dt, all_regions_anno, all.x = T, by = c("cluster_id"))
  return(dt)
}


#' Function to call DMRs on a camdac dmp dataset
#' @keywords internal
call_dmr_routine <- function(tmeth_dmps, regions_annotations, min_DMP_counts, min_consec_DMP) {
  # Annotate DMPs
  tmeth_dmps <- annotate_dmp_regions(tmeth_dmps, regions_annotations)
  dmp_cluster_counts <- get_cluster_counts(tmeth_dmps) # Note: Required later for merging to DMRs
  tmeth_dmps <- merge(tmeth_dmps, dmp_cluster_counts, by = "cluster_id", all.x = T)

  # Get CpG and DMP counts for each cluster
  tmeth_dmps[, DMR := ifelse(
    DMP_counts >= min_DMP_counts &
      consec_DMPs >= min_consec_DMP,
    "DMR",
    NA
  )]

  # Return NULL if no DMRs overlap
  # tmeth_dmr = tmeth_dmr[!is.na(DMR) & !is.na(DMP_t)]
  # if(nrow(tmeth_dmr) ==0){return(data.table())}

  # CAMDAC legacy: add CNA segment
  # tmeth_dmps[, segment := paste0(chrom,":",seg_start,"-",seg_end)]

  # Filter to DMRs
  tmeth_dmr <- collapse_cpg_to_dmr(tmeth_dmps)
  rm(tmeth_dmps)

  # Add CG counts and annotations
  tmeth_dmr <- merge(tmeth_dmr, dmp_cluster_counts, by = "cluster_id", all.x = T)
  tmeth_dmr <- re_annotate_dmrs(tmeth_dmr, regions_annotations)

  return(tmeth_dmr)
}



# Helper: Calculate the maximum number of consecutive DMPs in a dataset
max_consec_dmp <- function(x) {
  # Run-length encode DMPs
  rz <- rle(x[!is.na(x)])

  # Return NA if no DMPs present
  if (length(rz$lengths) == 0) {
    return(NA_integer_)
  }
  # Get maximum consecutive DMPs in DMR
  return(max(rz$lengths))
}

# Helper: Get DMR type
set_dmr_type <- function(x) {
  n_dmps <- length(x)
  is_hypo <- sum(x == "hypo") / n_dmps >= 0.9
  is_hyper <- sum(x == "hyper") / n_dmps >= 0.9
  is_mixed <- !is_hypo & !is_hyper
  dmr_type <- data.table::fcase(
    is_hypo, "hypo",
    is_hyper, "hyper",
    is_mixed, "mixed"
  )
  return(dmr_type)
}

call_dmrs <- function(tmeth_dmps, regions_annotations, itersplit = 3e5, min_DMP_counts = 5, min_consec_DMP = 4, n_cores = 5) {
  # Split region annotations in order to parallelise over subsets
  split_factor <- make_split_factor(nrow(regions_annotations), itersplit)
  regions_annotations <- split(regions_annotations, split_factor)

  #   Set warn=2 to ensure foreach fails if any of the parallel workers are terminated or raise a warning.
  #   without this option, foreach simply returns a warning and the pipeline continues. Essential for memory warning terminations.
  options(warn = 2)
  # Calculate DMR data for CpGs in parallel
  doParallel::registerDoParallel(cores = n_cores)
  dmrs <- foreach(regions_subset = regions_annotations, .combine = "rbind") %dopar% {
    call_dmr_routine(tmeth_dmps, regions_subset, min_DMP_counts, min_consec_DMP)
  }
  doParallel::stopImplicitCluster()
  options(warn = 0)

  return(dmrs)
}

# Note: mbdiff and mtdiff are calculated tumor - normal
dmp_call_pipe <- function(mbdiff, M_n, UM_n, M, UM, mtdiff = NULL, effect_size = 0.2, prob = 0.99, itersplit = 1e5, ncores=5) {
  # Calculate bulk DMP probability from counts
  phypo <- calc_prob_dmp(M_n, UM_n, M, UM, ncores = ncores, itersplit = itersplit)

  # Run bulk calculation if no pure given
  if (is.null(mtdiff)) {
    # If tumor is greater than normal then it's a hyper DMP
    prob_DMP <- data.table::fifelse(mbdiff > 0, 1 - phypo, phypo)
    DMP_b <- prob_to_call(prob_DMP, mbdiff, effect_size, prob)
    return(data.table(prob_DMP, mbdiff, DMP_b))
  } else {
    # Otherwise, run pure calculation and get DMP for bulk and pure
    prob_DMP <- data.table::fifelse(mtdiff > 0, 1 - phypo, phypo)
    DMP_b <- prob_to_call(prob_DMP, mbdiff, effect_size, prob)
    DMP_t <- prob_to_call(prob_DMP, mtdiff, effect_size, prob)
    return(data.table(prob_DMP, mbdiff, DMP_b, mtdiff, DMP_t))
  }
}
