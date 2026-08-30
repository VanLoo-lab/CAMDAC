# Main

cwrap_get_epialleles <- function(bam_file, seg, loci_dt = NA, paired_end, drop_ccgg,
                                 min_mapq = 1, min_cov = 3) {
    # Loci may be NA if loci for segment chromosome (i.e. chromY)
    # is missing. Return early in these cases as no alleles to count.
    # Two conditions required to avoid error raised using is.na alone.
    if (all(class(loci_dt) == "logical")) {
        if (is.na(loci_dt)) {
            return(empty_count_alleles_result())
        }
    }

    # Pre-applying multi-SNP loci filter
    loci_dt <- loci_dt[!duplicated(loci_dt, by = c("chrom", "start", "end"), fromLast = T)]

    # Read BAM and annotate SNP and CPG loci
    bam_dt <- get_reads_in_segments(bam_file, seg, min_mapq, paired_end = paired_end)
    if (nrow(bam_dt) == 0) {
        return(empty_count_alleles_result())
    }

    # Overlap with loci
    bam_dt <- format_bam_for_loci_overlap(bam_dt, paired_end = paired_end)
    bam_dt <- annotate_bam_with_loci(bam_dt, loci_dt, drop_ccgg = drop_ccgg, paired_end = paired_end)
    bam_dt <- drop_positions_outside_segments(bam_dt, seg)
    if (nrow(bam_dt) == 0) {
        return(empty_count_alleles_result())
    }

    if (paired_end) {
        # For paired end samples, we must select a single read at overlapping and then fix the
        # strand column to reflect read orientation (as per single-end CAMDAC)
        bam_dt <- fix_pe_overlap_at_loci(bam_dt) # Filters
        bam_dt <- add_loci_read_position(bam_dt)
        bam_dt <- fix_pe_strand_with_flags(bam_dt)
        bam_dt <- get_alleles_and_qual(bam_dt)
        bam_dt <- drop_pe_fields(bam_dt)
    } else {
        bam_dt <- add_loci_read_position(bam_dt)
        bam_dt <- get_alleles_and_qual(bam_dt)
    }

    # Additional filtering
    bam_dt <- filter_clipped_dinucleotides(bam_dt)
    # bam_dt <- filter_bam_by_quality(bam_dt, min_mapq = min_mapq)
    # Instead, annotate quality as we don't want to remove the knowledge of ref CG sites overlapping?
    hi_qual_dinucs <- data.table::like(bam_dt$qual.dinucs, "([5-9A-K:-@])([5-9A-K:-@])")
    hi_qual_snps <- data.table::like(bam_dt$qual.SNP, "([5-9A-K:-@])")

    # Set dinuc data to NA at positions where only SNP passed the quality filter
    #  This keeps the SNP data for downstream BAF/LogR but excludes the site from
    #  methylation rate calculations
    bam_dt[
        width >= 2 & !hi_qual_dinucs,
        c("alleles.dinucs", "qual.dinucs") := NA
    ]

    # Filter BAM for high quality dinucleotides and SNPs
    # bam_dt <- bam_dt[
    #   (width >= 2 & hi_qual_dinucs) | # Hi quality dinucleotide filter
    #     (!is.na(POS) & hi_qual_snps) # Hi quality SNP filter
    # ]

    # [1] 129852     20
    ## dim(bam_dt)

    # Filter records for minimum mapping quality
    # Note: mq filtering can also be applied by ScanBam
    bam_dt <- bam_dt[mq >= min_mapq]

    if (nrow(bam_dt) == 0) {
        return(empty_count_alleles_result())
    }

    # Get nucleotide counts and flatten pileup
    bam_dt <- annotate_nucleotide_counts(bam_dt)


    # RE filters: Up to this point, reads only filtered, not CpG sites. So all reads overlapping should be canonical.

    # At this stage, we have a mapping between reads and CpGs, with information on whether the site is a CG SNP
    # And the methylated/unmethylated counts
    # We don't care about SNPs on the read for now in general, so remove these
    bam_dt <- bam_dt[width > 1]

    # Now count the M or UM on the single CpGs, not on the pileup (as would be done in AC)
    bam_dt <- get_read_methylation_counts(bam_dt) # This function is defined here, above.

    bam_dt <- bam_dt[, chrom := gsub("chr", "", chrom)]

    # Apply CADMAC rules to get allele counts, methylation rates and BAFs
    res <- format_get_epialleles_result(bam_dt)

    return(res)
}

# Helper functions for entropy allele counter
format_get_epialleles_result <- function(x) {
    x[,
        .(
            read.start = min(read.start),
            read.end = min(read.end),
            cgstarts = list(start),
            states = list(rle(M))
        ),
        by = c("qname", "chrom")
    ]
}

get_read_methylation_counts <- function(bam_dt) {
    # A combination of get_snp_allele counts and get_read_methylation_counts
    # designed for single reads (not pileup)

    # Set REF_ALT combination string from loci
    bam_dt[!is.na(ref), SNP := paste0(ref, alt)]

    # Count REF alleles. Any missing loci are set to NA
    bam_dt[, ref_counts := data.table::fcase(
        SNP == "AC", Af + Ar,
        SNP == "CA", Tf + Cr + Cf,
        SNP == "AG", Af,
        SNP == "GA", Gf,
        SNP == "AT", Af + Ar,
        SNP == "TA", Tf + Tr,
        SNP == "GT", Gf + Ar + Gr,
        SNP == "TG", Tf + Tr,
        SNP == "CG", Tf + Cr + Cf,
        SNP == "GC", Gf + Ar + Gr,
        SNP == "CT", Cr,
        SNP == "TC", Tr
    )]

    # Count ALT alleles. Any loci not present are set to NA
    bam_dt[, alt_counts := data.table::fcase(
        SNP == "AC", Tf + Cr + Cf,
        SNP == "CA", Af + Ar,
        SNP == "AG", Gf,
        SNP == "GA", Af,
        SNP == "AT", Tf + Tr,
        SNP == "TA", Af + Ar,
        SNP == "GT", Tf + Tr,
        SNP == "TG", Gf + Ar + Gr,
        SNP == "CG", Gf + Ar + Gr,
        SNP == "GC", Tf + Cr + Cf,
        SNP == "CT", Tr,
        SNP == "TC", Cr
    )]
    # Count total reads contributing to SNP ref/alt counts
    bam_dt[, total_counts := ref_counts + alt_counts]

    # Count all reads with nucleotides expected by CAMDAC rules. This includes
    # positions where we couldn't distinguish between bisulfite conversion and SNPs.
    # This will later be subtracted from total depth to determine unexpected nucleotide count.
    bam_dt[, all_counts := data.table::fcase(
        is.na(ref), TGf + CAr + CGf + CGr,
        # For CT/AG SNPs, expected nucleotides are not in total_counts because they
        # confound bisulfite conversion, however we add them here for all_counts.
        SNP %like% "([GA][AG])", Af + Gf + Ar + Gr,
        SNP %like% "([CT][TC])", Cr + Tr + Tf + Cf,
        !is.na(ref), total_counts # All other positions get ref/alt counts
    )]

    # Set M from reads reporting methylation at CG dinucleotides,
    #  ignoring reads at CG-destroying SNPs.
    # Note: In RRBS version, CCGGs must be matched to SNP positions with a +1 offset
    bam_dt[
        width > 1, # Calculate for CpG loci only, setting non-CpGs to NA
        M := data.table::fcase(
            is.na(SNP), CGf + CGr,
            # When CpG starts at C/T SNP loci, we can't differentiate SNP from bisulfite conversion.
            # Therefore, count the reverse strand (bottom in directional library) only
            start == POS & SNP %like% "([CT][TC])", CGr,
            # When CpG ends at A/G SNP loci, we can't differentiate SNP from bisulfite conversion.
            # Therefore, count the forward strand (top in directional library) only
            end == POS & SNP %like% "([AG][GA])", CGf,
            # Count CpG dinucleotides at all other SNPs loci. This works as fcase
            # moves through conditions in order.
            !is.na(SNP), CGf + CGr
        )
    ]

    # Set UM as reads reporting unmethylated CG dinucleotides as above,
    # ignoring reads at CG-destroying SNPs.
    bam_dt[
        width > 1,
        UM := data.table::fcase(
            is.na(SNP), TGf + CAr,
            start == POS & SNP %like% "([CT][TC])", CAr,
            end == POS & SNP %like% "([AG][GA])", TGf,
            !is.na(SNP), TGf + CAr
        )
    ]

    # Set methylation counts to NA at loci with insufficient data
    bam_dt[
        M == 0 & UM == 0,
        c("M", "UM") := NA
    ]

    return(bam_dt)
}

filter_loci_dt_multi_snp <- function(x) {
    # Filter duplicated CpGs derived from CAMDAC annotations file
    x[, msl := duplicated(x, by = c("chrom", "start", "end"), fromLast = F) |
        duplicated(x, by = c("chrom", "start", "end"), fromLast = T)]
}

# Function to read and write epiallele states from a data.table
write_epi_states <- function(x, outfile, inplace = F) {
    if (inplace == F) { # Function below edits x in place. Decide whether to copy x or not
        x <- copy(x)
    }
    # x: data.table returned from entropy allele counter.
    # First, convert cg starts to pipe-separated string
    x[, cgstarts := lapply(cgstarts, function(x) paste0(x, collapse = "|"))]
    #  Next, convert states to pipe-separated string
    x[, states := lapply(states, function(x) paste0(inverse.rle(x), collapse = "|"))]
    # Finally, save to file
    data.table::fwrite(x, outfile)
}

read_epi_states <- function(infile) {
    x <- data.table::fread(infile)
    # First, convert cg starts to list
    x[, cgstarts := lapply(cgstarts, function(x) strsplit(x, "|", fixed = TRUE))]
    x[, cgstarts := lapply(cgstarts, function(x) as.numeric(x))]

    # Next, convert cgstates to list
    x[, states := lapply(states, function(x) strsplit(x, "|", fixed = TRUE))]
    suppressWarnings(x[, states := lapply(states, function(x) rle(as.numeric(x)))])

    return(x)
}

generate_epiallele_matrix <- function(x, chrom, region_start, region_end) {
    # x is an epiallele counts table
    # returns: A list with first entry as methylation matrix for the region, where
    # 0:unmethylated, 1:methylated, 2:unknown, and NA:unmapped
    # TODO: Do we want to restrict to reference genome hgs?
    # For example, CG-SNPs will always be NA and mess up entropy calculations
    stopifnot(region_end > region_start)

    x <- x[chrom == chrom]

    # Get reads that overlap region
    overlaps <- x[
        !(read.end < region_start | read.start > region_end)
    ]

    # Return empty list if no reads overlap
    if (nrow(overlaps) == 0) {
        return(list())
    }

    # Get the set of start sites in the matrix
    startsets <- unique(unlist(overlaps$cgstarts))
    st_order <- order(startsets)
    startsets <- startsets[st_order]
    ncgs <- length(startsets)

    # For each read, convert the states to a vector of length ncgs
    state_vecs <- lapply(
        seq(nrow(overlaps)),
        function(i) {
            # Get sites in the read's CpGs
            cur_rle <- overlaps$states[[i]]
            ssite <- inverse.rle(cur_rle)
            ssite[is.na(ssite)] <- 2 #  Sites uncovered are set to 2

            # Get vector of all sites in the region, setting to NA if not in read
            rsite <- overlaps$cgstarts[[i]]
            rsel <- match(rsite, startsets) # Vectorised which()
            cbool <- rep(NA, ncgs)
            cbool[rsel] <- ssite

            return(cbool)
        }
    )

    state_matrix <- do.call(rbind, state_vecs)
    colnames(state_matrix) <- startsets
    rownames(state_matrix) <- overlaps$qname

    # Keep CpGs in region
    col_sel <- startsets >= region_start & startsets <= region_end
    state_matrix <- state_matrix[, col_sel]

    res <- list(
        "matrix" = state_matrix,
        "meta" = data.frame(
            chrom = chrom,
            region_start = region_start,
            region_end = region_end,
            startsets = startsets
        )
    )

    return(res)
}
