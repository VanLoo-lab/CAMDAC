#' Count alleles for reads phased to SNPs in a BAM file
#' @param bam_file Path to BAM file
#' @param snps_gr GRanges object with heterozygous SNP loci for phasing
#' @param loci_dt Data table with CAMDAC CpG loci from reference files
#' @param paired_end Logical indicating if BAM is paired end
#' @param drop_ccgg Logical indicating if CCGG should be dropped (i.e. rrbs mode)
#' @param min_mapq Minimum mapping quality to consider a read
#' @param min_cov Minimum coverage to consider a read
#' @return A list with three slots: stats, qnames and asm_cg. stats describes counts of reads phased,
#'   qnames determines which SNPs each read was phased to and asm_cg is the data table with read counts
# Approx 1 minute for 5K SNPs on one core
cwrap_asm_get_allele_counts <- function(
    bam_file, snps_gr, loci_dt,
    paired_end, drop_ccgg, min_mapq = min_mapq, min_cov = min_cov) {

    # Read BAM
    bam_dt <- get_reads_in_segments(bam_file, snps_gr, min_mapq, paired_end = paired_end)
    bam_dt <- format_bam_for_loci_overlap(bam_dt, paired_end = paired_end)

    # Early strand adjustment for paired end: strands now reflect Watson/Crick strand (directional lib)
    bam_dt <- fix_pe_strand_with_flags(bam_dt, paired_end)

    # Overlap with SNP loci
    bam_dt <- phase_reads_to_snps(bam_dt, snps_gr)
    bam_dt <- select_read_snp_pair(bam_dt)

    # Assign alleles using CAMDAC rules
    bam_dt[, hap_is_ref := assign_het_allele(hap_bsseq, hap_ref, hap_alt, "ref")]
    bam_dt[, hap_is_alt := assign_het_allele(hap_bsseq, hap_ref, hap_alt, "alt")]

    # Get haplotype stats
    hap_stats <- asm_hap_stats(bam_dt)

    # Annotate BAM with loci
    bam_dt <- annotate_bam_with_loci_asm(bam_dt, loci_dt, drop_ccgg, paired_end)

    # Get qname to cpg mapping
    qname_hap_cg <- unique(bam_dt[, .(qname, hap_id, chrom, start, end)])

    # Count unexpected alleles and split by ref and alt after filtering
    ref_bam <- bam_dt[hap_is_ref == T]
    alt_bam <- bam_dt[hap_is_alt == T]
    rm(bam_dt)
    gc() # Remove bam and free up memory

    # Get counts for reads phased to each allele
    alt_cg <- asm_bam_to_counts(alt_bam, "alt", loci_dt,
        drop_ccgg = drop_ccgg, paired_end = paired_end
    )
    ref_cg <- asm_bam_to_counts(ref_bam, "ref", loci_dt,
        drop_ccgg = drop_ccgg, paired_end = paired_end
    )

    # Combine counts
    asm_cg <- merge(alt_cg, ref_cg,
        by = c(
            "CHR", "chrom", "start", "end",
            "width", "POS", "ref", "alt"
        ), all = TRUE
    )

    # Complete results object form hap_stats
    return(
        list(
            "asm_cg"=asm_cg,
            "hap_stats"=hap_stats,
            "map"=qname_hap_cg
        )
    )
    return(hap_stats)
}

haps_as_numeric <- function(v) {
    # v = c("1234;12", "1", "123;12;1")
    hap <- stringr::str_split(v, ";", simplify = T)
    hap <- as.numeric(hap)
    hap <- hap[!is.na(hap)]
    hap <- unique(hap)
    return(hap)
}

# FUTURE: Select based on counts
select_read_snp_pair <- function(bam_dt) {
    # Reads may map to multiple SNPs.
    # Ensure each read is represented only once.
    unique(bam_dt, by = "qname")
}

phase_reads_to_snps <- function(bam_dt, snps_gr) {
    # Find overlapping read pairs between BAM and SNPs
    snps_dt <- data.table(data.frame(snps_gr))
    setnames("seqnames", "chrom", x = snps_dt)
    setkey(bam_dt, chrom, start, end)
    setkey(snps_dt, chrom, start, end)
    bphase <- foverlaps(bam_dt, snps_dt, which = T, nomatch = NULL)
    setnames(bphase, "xid", "bam")
    setnames(bphase, "yid", "snps")

    # Run GAlignments parser on pairs to get the SNP position in reads
    bam_dt <- bam_dt[bphase$bam, ]
    snps_ph <- snps_dt[bphase$snps, ]
    aln <- GAlignments(
        seqnames = as.character(bam_dt$chrom), pos = bam_dt$start,
        cigar = as.character(bam_dt$cigar),
        strand = GenomicAlignments::strand(bam_dt$strand),
        names = as.character(bam_dt$qname)
    )
    gr <- GRanges(seqnames = snps_ph$chrom, ranges = IRanges(snps_ph$start, snps_ph$end))
    rpos <- pmapToAlignments(gr, aln)

    # Set haplotype information
    bam_dt$hap_ref <- snps_ph$ref
    bam_dt$hap_alt <- snps_ph$alt
    bam_dt$hap_id <- snps_ph$hap_id
    bam_dt$hap_POS <- snps_ph$start
    bam_dt$hap_allele <- substr(bam_dt$seq, start(rpos), end(rpos))
    bam_dt$hap_qual <- substr(bam_dt$qual, start(rpos), end(rpos))
    bam_dt$hap_bsseq <- paste0(bam_dt$hap_allele, bam_dt$strand)
    return(bam_dt)
}

assign_het_allele <- function(bseq_strand, ref, alt, call) {
    stopifnot(call %in% c("ref", "alt"))
    SNP <- paste0(ref, alt)
    if (call == "ref") {
        bool <- dplyr::case_when(
            SNP == "AC" & (bseq_strand %in% c("A+", "A-")) ~ TRUE,
            SNP == "CA" & (bseq_strand %in% c("T+", "C-", "C+")) ~ TRUE,
            SNP == "AG" & (bseq_strand %in% c("A+")) ~ TRUE,
            SNP == "GA" & (bseq_strand %in% c("G+")) ~ TRUE,
            SNP == "AT" & (bseq_strand %in% c("A+", "A-")) ~ TRUE,
            SNP == "TA" & (bseq_strand %in% c("T+", "T-")) ~ TRUE,
            SNP == "GT" & (bseq_strand %in% c("G+", "A-", "A+")) ~ TRUE,
            SNP == "TG" & (bseq_strand %in% c("T+", "T-")) ~ TRUE,
            SNP == "CG" & (bseq_strand %in% c("T+", "C-", "C+")) ~ TRUE,
            SNP == "GC" & (bseq_strand %in% c("G+", "A-", "G-")) ~ TRUE,
            SNP == "CT" & (bseq_strand %in% c("C-")) ~ TRUE,
            SNP == "TC" & (bseq_strand %in% c("T-")) ~ TRUE,
            TRUE ~ FALSE
        )
    }


    if (call == "alt") {
        bool <- dplyr::case_when(
            SNP == "AC" & (bseq_strand %in% c("T+", "C-", "C+")) ~ TRUE,
            SNP == "CA" & (bseq_strand %in% c("A+", "A-")) ~ TRUE,
            SNP == "AG" & (bseq_strand %in% c("G+")) ~ TRUE,
            SNP == "GA" & (bseq_strand %in% c("A+")) ~ TRUE,
            SNP == "AT" & (bseq_strand %in% c("T+", "T-")) ~ TRUE,
            SNP == "TA" & (bseq_strand %in% c("A+", "A-")) ~ TRUE,
            SNP == "GT" & (bseq_strand %in% c("T+", "T-")) ~ TRUE,
            SNP == "TG" & (bseq_strand %in% c("G+", "A-", "G-")) ~ TRUE,
            SNP == "CG" & (bseq_strand %in% c("G+", "A-", "G-")) ~ TRUE,
            SNP == "GC" & (bseq_strand %in% c("C-", "T+", "C+")) ~ TRUE,
            SNP == "CT" & (bseq_strand %in% c("T-")) ~ TRUE,
            SNP == "TC" & (bseq_strand %in% c("C-")) ~ TRUE,
            TRUE ~ FALSE
        )
    }
    return(bool)
}

# Summarise read counts on haplotypes after CAMDAC allele counting rules
asm_hap_stats <- function(bam_dt) {
    # Get count of sites that could not be included in counts
    # These represent unexpected nucleotides (e.g. SNV) and sites where bisulfite leaves ambiguous
    bam_dt[hap_is_ref == F & hap_is_alt == F, hap_unexp := 1]
    unexp_dt <- bam_dt[, .(hap_unexp = sum(hap_unexp, na.rm = T)), by = c("chrom", "hap_ref", "hap_alt", "hap_POS", "hap_id")]

    # Select reads that would be taken for downstream analysis
    bam_dt <- bam_dt[hap_is_ref == T | hap_is_alt == T]

    # Count reads aligned to input haplotype/SNP
    stats <- bam_dt[, .(hap_BAF = sum(hap_allele == hap_alt) / .N, hap_reads = .N), by = c("chrom", "hap_ref", "hap_alt", "hap_POS", "hap_id")]
    stats <- merge(stats, unexp_dt, all.x = T)

    # Ensure BAM chrom field fits expected format for downstream joins
    stats$chrom <- gsub("chr", "", stats$chrom)

    # Return stats
    return(stats)
}

asm_bam_to_counts <- function(
    asm_dt, asm_type, loci_dt, drop_ccgg = FALSE,
    paired_end = FALSE, min_mapq = 0) {
    stopifnot(asm_type %in% c("ref", "alt"))

    if (paired_end) {
        asm_dt <- fix_pe_overlap_at_loci(asm_dt)
        # N.B. pe strand fixed earlier in pipeline
        asm_dt <- add_loci_read_position(asm_dt)
        asm_dt <- get_alleles_and_qual(asm_dt)
        asm_dt <- drop_pe_fields(asm_dt)
    } else {
        asm_dt <- add_loci_read_position(asm_dt)
        asm_dt <- get_alleles_and_qual(asm_dt)
    }

    # Additional filtering
    asm_dt <- filter_clipped_dinucleotides(asm_dt)
    asm_dt <- filter_bam_by_quality(asm_dt, min_mapq = min_mapq)

    # Pileup
    asm_dt <- annotate_nucleotide_counts(asm_dt)
    pileup_summary <- flatten_pileup_to_counts(asm_dt)
    rm(asm_dt)

    # Apply CADMAC rules to get allele counts, methylation rates and BAFs
    pileup_summary <- get_snp_allele_counts(pileup_summary)
    pileup_summary <- get_methylation_counts(pileup_summary, min_cov)
    pileup_summary <- filter_bad_allele_count_rows(pileup_summary, min_cov)
    pileup_summary <- compute_methylation_rates(pileup_summary)
    pileup_summary <- compute_BAFs(pileup_summary)

    # Format result for output
    result <- format_get_reads_result(pileup_summary)

    # Drop any sites with SNPs counts only but no methylation
    result <- result[!is.na(total_counts_m), ]

    # Select fields to keep
    result <- result[, .(
        CHR, chrom, start, end, width, POS, ref, alt,
        alt_counts, ref_counts, total_counts, BAF, total_depth,
        other_counts, all_counts, M, UM, total_counts_m, m, CCGG
    )]

    rename_cols <- setdiff(
        names(result),
        c(
            "CHR", "chrom", "start", "end", "width", "POS", "ref",
            "alt"
        )
    )
    # Give ref/alt names to essential columns
    for (n in rename_cols) {
        setnames(result, n, paste0(asm_type, "_", n))
    }
    return(result)
}

write_asm_counts_output <- function(result, sample, config){
  cg_outfile <- get_fpath(sample, config, "asm_counts")
  data.table::fwrite(result$asm_cg, outfile)

  phase_outfile <- get_fpath(sample, config, "asm_phase_map")
  data.table::fwrite(result$map, phase_outfile) 
  
  stats_outfile <- get_fpath(sample, config, "asm_hap_stats")
  data.table::fwrite(result$hap_stats, stats_outfile)

  return(cg_outfile)
}

annotate_bam_with_loci_asm <- function(bam_dt, loci_subset, drop_ccgg=F, paired_end=F){
    # Set keys for join
    loci_subset$chrom = as.character(loci_subset$chrom)
    data.table::setkey(loci_subset, chrom, start, end)
    bam_dt$chrom = as.character(bam_dt$chrom)
    data.table::setkey(bam_dt, chrom, start, end)

    # Filter CCGG loci if WGBS
    if (drop_ccgg) {
    loci_subset <- loci_subset[width != 4]
    }

    # Overlap
    bam_loci_overlap <- data.table::foverlaps(bam_dt, loci_subset)

    # Rename read fields
    setnames(bam_loci_overlap, "i.start", "read.start")
    setnames(bam_loci_overlap, "i.end", "read.end")
    setnames(bam_loci_overlap, "mapq", "mq")
    bam_loci_overlap[, strand := i.strand] # Set strand as previous

    # Filter out rows with no loci data, set expected columns and return
    bam_loci_overlap <- bam_loci_overlap[!is.na(width), ]

    return(bam_loci_overlap)
}

load_asm_loci_for_segment <- function(snps_gr, loci_files){
    snps_region = reduce(snps_gr+1000) # Get regions in 1kb non-overlapping regions around SNPs
    loci_dt <- load_loci_for_segment(snps_region, loci_files)
    loci_dt <- loci_dt[width > 1, ] # Ensure only CG sites are mapped for ASM 
    return(loci_dt)
}