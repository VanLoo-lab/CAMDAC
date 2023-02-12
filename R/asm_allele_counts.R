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
cwrap_asm_get_allele_counts <- function(bam_file, snps_gr, loci_dt,
paired_end, drop_ccgg, min_mapq = min_mapq, min_cov = min_cov) {
    # Ensure only CGs
    loci_dt <- loci_dt[width > 1, ]

    # Read BAM
    bam_dt <- get_reads_in_segments(bam_file, snps_gr, min_mapq, paired_end = paired_end)
    bam_dt <- format_bam_for_loci_overlap(bam_dt, paired_end = paired_end)

    # Overlap with SNP loci
    bam_dt <- phase_reads_to_snps(bam_dt, snps_gr)
    bam_dt <- select_read_snp_pair(bam_dt)

    # Early strand adjustment for paired end: strands now reflect Watson/Crick strand (directional lib)
    bam_dt <- fix_pe_strand_with_flags(bam_dt, paired_end)

    # Assign alleles using CAMDAC rules
    bam_dt[, hap_is_ref := assign_het_allele(hap_bsseq, hap_ref, hap_alt, "ref")]
    bam_dt[, hap_is_alt := assign_het_allele(hap_bsseq, hap_ref, hap_alt, "alt")]

    # Get haplotype stats for downstream analysis
    hap_stats <- asm_hap_stats(bam_dt)

    # Count unexpected alleles and split by ref and alt after filtering
    ref_bam <- bam_dt[hap_is_ref == T]
    alt_bam <- bam_dt[hap_is_alt == T]
    rm(bam_dt);gc(); # Remove bam and free up memory

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


    # Get haplotype data as numeric
    hap_ix <- list(
        "ref"=haps_as_numeric(asm_cg$ref_hap_cg),
        "alt"=haps_as_numeric(asm_cg$alt_hap_cg)
    )

    # Complete results object form hap_stats
    hap_stats$asm_cg = asm_cg
    hap_stats$hap_ix = hap_ix

    return(hap_stats)
}

haps_as_numeric <- function(v){
    # v = c("1234;12", "1", "123;12;1")
    hap = stringr::str_split(v, ';', simplify=T)
    hap = as.numeric(hap)
    hap = hap[!is.na(hap)]
    return(hap)
}

# FUTURE: Select based on counts
select_read_snp_pair <- function(bam_dt){
    # Reads may map to multiple SNPs.
    # Ensure each read is represented only once.
    unique(bam_dt, by="qname")
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
    gr <- GRanges(seqnames = snps_ph$chrom, ranges = IRanges(snps_ph$POS, snps_ph$POS))
    rpos <- pmapToAlignments(gr, aln)

    # Set haplotype information
    bam_dt$hap_ref <- snps_ph$ref
    bam_dt$hap_alt <- snps_ph$alt
    bam_dt$hap_id <- snps_ph$hap_id
    bam_dt$hap_POS <- snps_ph$POS
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
    bam_dt[hap_is_ref == F & hap_is_alt == F, hap_unexp:=1]
    unexp_dt = bam_dt[, .(hap_unexp = sum(hap_unexp, na.rm=T)), by = c("chrom", "hap_ref", "hap_alt", "hap_POS", "hap_id")]

    # Select reads that would be taken for downstream analysis
    bam_dt <- bam_dt[hap_is_ref == T | hap_is_alt == T]

    # Count reads aligned to input haplotype/SNP
    stats <- bam_dt[, .(hap_BAF = sum(hap_allele == hap_alt) / .N, hap_reads = .N), by = c("chrom", "hap_ref", "hap_alt", "hap_POS", "hap_id")]
    stats <- merge(stats, unexp_dt, all.x=T)

    # Ensure BAM chrom field fits expected format for downstream joins
    stats$chrom <- gsub("chr", "", stats$chrom)

    # Get selection of qname to hap_id mapping (1 to 1 for each read)
    qname_hap_id <- unique(bam_dt[, .(qname, hap_id)])

    # Return stats
    obj <- list(stats=stats, qnames=qname_hap_id)
    return(obj)
}


empty_asm_bam_to_counts <- function(asm_type){
     # Set empty return value
    fields <- c("CHR", "chrom", "start", "end", "width", "POS", 
    "ref", "alt")
    count_names <- c(
        "alt_counts", "ref_counts",
        "total_counts", "BAF", "total_depth", "other_counts", "all_counts",
        "M", "UM", "total_counts_m", "m", "CCGG"
    )
    empty_return <- data.table(matrix(nrow=0, ncol=length(c(fields, count_names))))
    names(empty_return) = c(fields, paste0(asm_type, "_", count_names))

    # Set field classes for downstream joins
    empty_return$chrom = as.character(empty_return$chrom)
    empty_return$ref = as.character(empty_return$ref)
    empty_return$alt = as.character(empty_return$alt)

    empty_return
}

asm_bam_to_counts <- function(
    asm_dt, asm_type, loci_dt, drop_ccgg = FALSE,
    paired_end = FALSE, min_mapq = 0) {
    stopifnot(asm_type %in% c("ref", "alt"))

    # Fix hap_id for downstream overlap
    hap_id_data = asm_dt[, .(qname, hap_id)]

    # Annotate BAM with CpG-SNP loci
    # As this is the main camdac allele counter annotator,
    # the hap_id field is removed by default
    asm_dt <- annotate_bam_with_loci(asm_dt, loci_dt,
        drop_ccgg = drop_ccgg, paired_end = paired_end
    )

    # Combine the BAM:CpG mappings with the BAM:hap_id map
    # CpGs may have reads mapped to multiple SNPs,
    # so combine their IDs with ';'
    hap_id_cg <- merge(
        asm_dt[,.(qname, chrom, start, end)],
        hap_id_data,
        all.x=T)[, 
        .(hap_cg=paste0(unique(hap_id), collapse=";")),
        by=c("chrom", "start", "end")
        ]
    hap_id_cg$chrom = gsub("chr", "", hap_id_cg$chrom)

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

    # Add hap_id so that we can join to haplotype stats
    result <- merge(result, hap_id_cg, by=c("chrom","start","end"), all.x=T)

    rename_cols = setdiff(
        names(result),
        c("CHR", "chrom", "start", "end", "width", "POS", "ref",
        "alt")
    )
    # Give ref/alt names to essential columns
    for (n in rename_cols) {
        setnames(result, n, paste0(asm_type, "_", n))
    }
    return(result)
}
