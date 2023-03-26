overlap_meth_cna <- function(asm_meth, cna, hap_stats, phase_map) {
    # Format data for merge
    setkey(asm_meth, chrom, start, end)
    setkey(cna, chrom, start, end)

    # Annotate BAF of SNP for each phased CpG
    baf <- hap_stats[, .(
        hap_POS, hap_id, hap_BAF, hap_reads,
        ref_allele = hap_ref, alt_allele = hap_alt
    )]
    phase <- unique(phase_map[, .(chrom, start, end, hap_id)])
    phase$chrom <- stringr::str_replace(phase$chrom, "chr", "")
    phase <- merge(phase, baf, all.y = T, by = "hap_id") # all.y required?

    # Merge methylation and cna
    amc <- data.table::foverlaps(asm_meth, cna, mult = "first", nomatch = NA)
    amc$start <- amc$cna_start
    amc$end <- amc$cna_end
    amc$start <- amc$i.start
    amc$end <- amc$i.end # Set back
    amc$i.start <- NULL
    amc$i.end <- NULL
    setkey(amc, chrom, start, end)

    # Merge meth_cna and BAF of phased SNPs, along with cg-phase annotation
    setkey(phase, chrom, start, end)
    res <- merge(amc, phase, all.x = T, by = c("chrom", "start", "end"))

    return(res)
}

assign_asm_cna <- function(ol) {
    # TODO: This is a naive and incorrect solution. Update assignment

    ol[
        nA == nB,
        `:=`(
            ref_CN = nA,
            alt_CN = nB
        )
    ]

    ol[ hap_BAF <= 0.5, `:=`(
        ref_CN = nA,
        alt_CN = nB
    )]

    ol[ hap_BAF > 0.5, `:=`(
        ref_CN = nB,
        alt_CN = nA
    )]
}
