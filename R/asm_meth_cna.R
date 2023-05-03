merge_asm_hap <- function(asm_meth, hap_stats, phase_map) {
    # Annotate BAF of SNP for each phased CpG
    baf <- hap_stats[, .(
        hap_POS, hap_id, hap_BAF, hap_reads,
        ref_allele = hap_ref, alt_allele = hap_alt
    )]
    phase <- unique(phase_map[, .(chrom, start, end, hap_id)])
    phase$chrom <- stringr::str_replace(phase$chrom, "chr", "")
    phase <- merge(phase, baf, all.y = T, by = "hap_id")

    # Merge methylation and haplotype stats
    setkey(phase, chrom, start, end)
    res <- merge(asm_meth, phase, all.x = T, by = c("chrom", "start", "end"))
    return(res)
}

overlap_meth_cna <- function(asm_hap, cna) {
    # Format data for merge
    setkey(asm_hap, chrom, start, end)
    setkey(cna, chrom, start, end)

    # Merge methylation and cna
    amc <- data.table::foverlaps(asm_hap, cna, mult = "first", nomatch = NA)
    amc$cna_start <- amc$start
    amc$cna_end <- amc$end
    amc$start <- amc$i.start # Set back to CG
    amc$end <- amc$i.end # Set back to CG
    amc$i.start <- NULL
    amc$i.end <- NULL
    setkey(amc, chrom, start, end)

    return(amc)
}

assign_asm_cna <- function(ol) {
    # TODO: Use battenberg phasing where available
    # TODO: Explain CpGs where ref_CN > alt_CN but ref_m is NA?
    ol = data.table(ol)

    # Simple case: equal CNA states
    ol[
        major_cn == minor_cn,
        `:=`(
            ref_CN = major_cn,
            alt_CN = minor_cn
        )
    ]

    # For each segment, annotate whether nMajor is consistent with ref or alt.
    segmap = ol[, .(
        mBAF = mean(hap_BAF, na.rm=TRUE),
        medBaF = median(hap_BAF, na.rm=TRUE),
        MajRef = sum(hap_BAF <= 0.5, na.rm=TRUE),
        MajAlt = sum(hap_BAF >= 0.5, na.rm=TRUE),
        nBAF = .N),
        by = c("chrom", "cna_start", "cna_end", "major_cn", "minor_cn")
    ]

    # Overlap
    ol <- left_join(ol, segmap)


    # Assign based on BAF per segment
    ol[
        is.na(ref_CN) & is.na(alt_CN) & (MajRef > MajAlt),
        `:=`(
            ref_CN = major_cn,
            alt_CN = minor_cn
        )
    ]
    ol[
        is.na(ref_CN) & is.na(alt_CN) & (MajRef < MajAlt),
        `:=`(
            ref_CN = minor_cn,
            alt_CN = major_cn
        )
    ]

    return(ol)
}
