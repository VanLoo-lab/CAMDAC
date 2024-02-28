# Import calls make package functions available to CAMDAC internal functions without :: syntax.
# NULL required for roxygen to generate these import calls in NAMESPACE after devtools::document()

#' @import data.table
#' @import ggplot2
#' @import foreach
#' @import doParallel
#' @import GenomicRanges
#' @import GenomicAlignments
#' @import logging
#' @import ggplot2
#' @import dplyr
#' @import gridExtra
#' @import stringr
#' @import png
#' @import MASS
#' @import GenomeInfoDb
#' @importFrom IRanges IRanges
#' @import S4Vectors
#' @import Rsamtools
#' @importFrom grDevices adjustcolor dev.off png rgb
#' @importFrom graphics abline axis par plot points rect text
#' @importFrom stats cor frequency lm median na.omit optimize qbeta rbeta runif setNames
#' @importFrom utils read.table write.table
NULL

# Global variables used in packages like data.table with non-standard evaluation.
# utils::globalVariables stops errors in devtools::check() and submissions to CRAN.
utils::globalVariables(
  c(
    "ref.snp",
    "ref",
    "alt.snp",
    "alt.af",
    "chrom",
    "start",
    "end",
    "nA",
    "nB",
    "seg",
    "i",
    "DMR",
    "DMP_t",
    "m_n",
    "m_x_low_n",
    "m_x_high_n",
    "m_t",
    "m_t_low",
    "m_t_high",
    "prob_DMP",
    "m",
    "CG_CN",
    "nA",
    "nB",
    "segment",
    "chrom",
    "i.start",
    "i.end",
    "ref",
    "alt",
    "total_counts",
    "alt_counts",
    "LogR_corr",
    "LogR",
    "GC",
    "repli",
    "POS",
    "ref_counts",
    "dt_chrom",
    "#CHR",
    "Count_A",
    "Count_C",
    "Count_G",
    "Count_T",
    "Good_depth",
    "..outfile_columns",
    "normalLogR",
    "p",
    "cov",
    "cov_t",
    "m_t",
    "end",
    "width",
    "read.start",
    "CN",
    "m_t_raw",
    "groupid",
    "N",
    "flag",
    "Af", "Ar", "Cf", "Cr", "Gf", "Gr", "Tf", "Tr", "CGf", "CGr", "TGf", "CAr", "CCGG", "mq",
    "POS", "qwidth", "..keep_columns", "CHR", "alleles.dinucs", "qual.dinucs", "gr_normal",
    "gr_tumor", "M", "UM", "total_counts_m", "M_snuc", "mate_status", "total_depth",
    "BAFloci", ".", "flag", "Frequency", "BAF", "cov", "rstart", "rend", "rwidth",
    "snp_width", "qname", "ah_strand", "snp_start", "strand", "cigar", "ID", "QUAL",
    "FILTER", "INFO", "FORMAT", "strand", "i.strand", "..bam_cols", "cg_snp",
    "CG_CN_n", "cluster_id", "seqnames", "gc_file", "BAFr", "total_counts_n", "alleles.SNP",
    "..cell_line_cols", "AF", "qual", "qual.SNP", "SNP", "UM_snuc", "all_counts", "other_counts",
    "selector", "nMaj1_A", "nMaj2_A", "nMin1_A", "nMin2_A", "chr", "startpos", "endpos",
    "is.best", "..bb_cna_fields", "rBAF", "BAF_n", "..count.."
  )
)