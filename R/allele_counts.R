detect_genome_build <- function(bam_file) {
  # Get genome build from BAM header for opening reference files

  # Get reference sequence lengths from BAM header
  ref_lengths <- Rsamtools::scanBamHeader(bam_file)[[1]][[1]]
  # Test whether chromosome X length is expected for hg19
  # X used as can be selected regardless of reference
  x_is_hg19 <- ref_lengths[grepl("X$", names(ref_lengths))] == "155270560"
  # Set build
  build <- ifelse(unname(x_is_hg19), "hg19", "hg38")

  return(build)
}

get_reads_in_segments <- function(bam_file, segments, min_mapq, paired_end = FALSE) {
  # Get reads from pre-determined regions (segments)

  # Set BAM flags required for paired end reads
  if (paired_end) {
    # For paired-end sequencing, we use the settings below to limit to following flags to ensure we
    # are only left with Bismark informative paired reads:
    #      R1  R2
    #  OT  99  147
    #  OB  83  163
    flag <- Rsamtools::scanBamFlag(
      isPaired = TRUE, isProperPair = TRUE,
      isSecondaryAlignment = FALSE, isSupplementaryAlignment = FALSE,
      isUnmappedQuery = FALSE, isDuplicate = FALSE
    )
    # Adding flag to determine PE OT/OB downstream
    what <- c("qname", "rname", "flag", "strand", "pos", "qwidth", "mapq", "seq", "qual", "groupid", "cigar", "mate_status")
    asMates <- TRUE # Adds additional mates column to for filtering. BamFile object only.
  } else {
    flag <- Rsamtools::scanBamFlag()
    what <- c("qname", "rname", "strand", "pos", "qwidth", "mapq", "seq", "qual", "cigar")
    asMates <- FALSE
  }

  # Set parameters and read BamFile.
  #    ScanBamParam "whats" options are available via scanBamWhat()
  param <- Rsamtools::ScanBamParam(which = segments, what = what, flag = flag, mapqFilter = min_mapq)
  bam_object <- Rsamtools::BamFile(bam_file, index = paste0(bam_file, ".bai"), asMates = asMates)
  # Read BAM file as a list of lists
  bam <- Rsamtools::scanBam(bam_object, param = param)

  # Convert BAM to data.table with one row per read
  # For standard WGBS analysis: 
  # Segments of the cmain function is generated from reference files or bed regions
  # Bam will have length of one as only one range is provided
  # For ASM, segments = snps_gr is a GRanges object with multiple ranges
  # Bam will have length > 1 which needs rbindlist and as.data.frame.list

  # The lapply operation takes bam from a list-of-lists to a list-of-dataframes,
  # which rbindlist combines into a single datatable.
  bam_dt <- data.table::rbindlist(
    lapply(bam, as.data.frame.list)
  )

  # UCSC Fix (Vectorized is faster)
  if (nrow(bam_dt) > 0 && !startsWith(as.character(bam_dt$rname[1]), "chr")) {
    bam_dt[, rname := paste0("chr", rname)]
  }

  return(bam_dt)

}

format_bam_for_loci_overlap <- function(bam_dt, paired_end = FALSE) {
  setnames(bam_dt, "rname", "chrom")
  setnames(bam_dt, "pos", "start")
  bam_dt[, end := (start + (qwidth - 1))]

  # Ensure chrom column is UCSC format
  #  CAMDAC-RRBS performs Ensembl/UCSC format checks on seqnames.
  #  CAMDAC-WGBS will use UCSC as default (e.g. chr1) except for ASCAT where Ensembl format is required
  chr_entry <- as.character(bam_dt[["chrom"]][[1]])
  if (!startsWith(chr_entry, "chr")) {
    bam_dt[, chrom := paste0("chr", chrom)]
  }

  keep_columns <- c("qname", "chrom", "strand", "start", "end", "mapq", "seq", "qual", "cigar")

  if (paired_end) {
    keep_columns <- c(keep_columns, "flag", "groupid", "mate_status")
  }

  # Return data table with columns filtered
  return(bam_dt[, ..keep_columns])
}

load_loci_as_data_table <- function(loci_file, drop_ccgg = TRUE) {
  base::load(loci_file) # Brings loci_subset into environment
  loci_dt <- data.table::data.table(data.frame(loci_subset))
  # Filter any CCGG sites out for WGBS analysis
  if (drop_ccgg) {
    loci_dt <- loci_dt[width != 4]
  }
  setnames(loci_dt, "seqnames", "chrom")
  return(loci_dt)
}

annotate_bam_with_loci <- function(bam_dt, loci_subset, paired_end = FALSE) {
  # DO NOT modify loci_subset here! Avoid deep copy in workers.
  # Ensure bam_dt matches the key structure of the pre-keyed loci
  bam_dt[, chrom := as.character(chrom)]
  data.table::setkey(bam_dt, chrom, start, end)
  bam_loci_overlap <- data.table::foverlaps(bam_dt, loci_subset, nomatch = NULL)

  # Filter overlap to expected columns and rename columns to those used in CAMDAC-RRBS
  data.table::setnames(bam_loci_overlap, 
                     old = c("i.start", "i.end", "mapq"), 
                     new = c("read.start", "read.end", "mq"))
  # Replace the i.strand column with the strand column
  # The strand column is "*" as it derives from loci_subset,
  # while i.strand derives from the BAM and contains true strand orientation for each read
  bam_loci_overlap[, `:=`(
    strand = i.strand,
    i.strand = NULL
  )]
  bam_cols <- c(
    "qname", "strand", "chrom", "read.start", "read.end", "POS", "width",
    "start", "end", "ref", "alt", "seq", "qual", "mq", "cigar"
  )
  if (paired_end) {
    bam_cols <- c(bam_cols, "flag", "groupid", "mate_status")
  }

  # Set expected columns and return
  return(bam_loci_overlap[, ..bam_cols])
}

drop_positions_outside_segments <- function(bam_dt, segments) {
  # Assuming segments will be an arbitrary number of ranges for a single chrom,
  # bind segments in gr with range() and get the start and end
  segs_range <- range(segments)
  segs_start <- start(segs_range)
  segs_end <- end(segs_range)

  # Drop bam_dt positions outside of these ranges
  bam_dt <- bam_dt[(start >= segs_start) & (start <= segs_end)]

  return(bam_dt)
}

infer_ot_ob_strand_from_flag <- function(flag) {
  # Check for OT/OB patterns using bitwise logic
  # Bismark OT/OB flags:
  #   OT R1: 64 (R1) + 32 (mate_rev) + 2 (proper) + 1 (paired) = 99
  #   OT R2: 128 (R2) + 16 (self_rev) + 2 + 1 = 147
  #   OB R1: 64 (R1) + 16 (self_rev) + 2 + 1 = 83
  #   OB R2: 128 (R2) + 32 (mate_rev) + 2 + 1 = 163

  # Bit definitions (SAM spec)
  is_paired      <- bitwAnd(flag, 0x1) != 0     # 1
  is_proper_pair <- bitwAnd(flag, 0x2) != 0     # 2
  is_reverse     <- bitwAnd(flag, 0x10) != 0    # 16
  mate_reverse   <- bitwAnd(flag, 0x20) != 0    # 32
  is_R1          <- bitwAnd(flag, 0x40) != 0    # 64
  is_R2          <- bitwAnd(flag, 0x80) != 0    # 128

  is_OT <- (is_R1 & !is_reverse & mate_reverse & is_proper_pair & is_paired) |
           (is_R2 & is_reverse & !mate_reverse & is_proper_pair & is_paired)

  is_OB <- (is_R1 & is_reverse & !mate_reverse & is_proper_pair & is_paired) |
           (is_R2 & !is_reverse & mate_reverse & is_proper_pair & is_paired)

  return(
    ifelse(is_OT, "+",
           ifelse(is_OB, "-", NA_character_))
  )
}

fix_pe_strand_with_flags <- function(bam_dt, paired_end = TRUE) {
  if (paired_end) {
  # Convert "strand" column to CAMDAC-expected strand using Bismark flags for OT/OB
  #
  #                 |  R1     R2    |  CAMDAC strand column
  # Bismark flag OT | 99(+)  147(-) |       OT = "+"
  # Bismark flag OB | 83(-)  163(+) |       OB = "-"
  # Note, this is the same as viewing in IGV as "first of pair strand".

    # DO NOT setkey(bam_dt, groupid). Sorting 20M rows here is a waste of time and RAM.
    # Use fcase (Fast Case) to map the specific Bismark OT/OB flags directly
    # 99/147 are CAMDAC Top (+), 83/163 are CAMDAC Bottom (-)
    bam_dt[, strand := data.table::fcase(
      flag == 99L  | flag == 147L, "+",
      flag == 83L  | flag == 163L, "-",
      default = NA_character_
    )]
  }
  return(bam_dt)
}

fix_pe_overlap_at_loci <- function(bam_dt) {
    # Filter mate pairs that overlap loci.
  # We prefentially keep R2 as for Swift Accel-MethylSeq library, the tail of +ve R1 may
  # contain adapter contaminant sequences typically trimmed off 5' +ve R2.
  # See: https://swiftbiosci.com/wp-content/uploads/2019/02/16-0853-Tail-Trim-Final-442019.pdf
  # 210401 - I found that the R1/"first-in-pair" flags actually have better per-base quality
  #    than their R2 counterparts. As we aren't using PE overlaps, filter R2 to see if it improves score

  # Set flag for R1 status
  data.table::setkey(bam_dt, NULL)
  is_r1 <- bitwAnd(bam_dt$flag, 0x40) != 0
  
  dup_cols <- c("chrom", "width", "start", "groupid")
  is_dup <- duplicated(bam_dt, by = dup_cols) | 
            duplicated(bam_dt, by = dup_cols, fromLast = TRUE)

  bam_dt <- bam_dt[!is_dup | is_r1]
  return(bam_dt)
}

add_loci_read_position_skipCIGAR <- function(bam_dt) {
  # Take only reads without indels or clipping
  bam_dt <- bam_dt[grepl("^\\d+M$", cigar)]

  # Set rstart and rend by parsing relative to genome
  # This depends on strand orientation simply because read start and end are reversed
  bam_dt[, `:=`(
    rstart = (start - read.start + 1),
    rend = (end - read.start + 1)
  )]

  # Set ccgg sites to 0 for now
  # TODO: Process PE CCGG (RRBS) in CAMDAC?
  bam_dt[width == 4, `:=`(rstart = 0, rend = 0)]

  return(bam_dt)
}

add_loci_read_position <- function(bam_dt) {
  # 1. Identify 'Simple' reads: Only digits followed by 'M' (e.g., 150M)
  # This avoids calling the heavy CIGAR parser for 99% of the data
  # We use fixed=FALSE but a very specific regex for speed.
  is_complex <- grepl("[IDSHNX=]", bam_dt$cigar, perl = TRUE)
  
  # 2. For simple 'M' matches, the read position is just the genomic 
  # offset from the start of the alignment.
  bam_dt[is_complex == FALSE, `:=`(
    rstart = as.integer(start - read.start + 1L),
    rend = as.integer(end - read.start + 1L)
  )]
  
  # 3. Only for the reads with Indels/Clipping
  if (any(is_complex)) {
    complex_idx <- which(is_complex)
    sub_dt <- bam_dt[complex_idx]
    
    aln <- GenomicAlignments::GAlignments(
      seqnames = as.character(sub_dt$chrom), 
      pos = sub_dt$read.start,
      cigar = as.character(sub_dt$cigar), 
      strand = GenomicAlignments::strand(sub_dt$strand),
      names = as.character(sub_dt$qname)
    )
    
    gr <- GenomicRanges::GRanges(
      seqnames = sub_dt$chrom, 
      ranges = IRanges::IRanges(sub_dt$start, sub_dt$end)
    )
    
    res <- GenomicAlignments::pmapToAlignments(gr, aln)
    bam_dt[complex_idx, `:=`(
      rstart = GenomicRanges::start(res), 
      rend = GenomicRanges::end(res)
    )]
  }

  keep_idx <- bam_dt[, which(
    rstart >= 1L & 
    rend <= as.integer(read.end - read.start + 1L) & 
    width == (rend - rstart + 1L)
  )]
  
  return(bam_dt[keep_idx])
}

add_loci_read_position_legacy <- function(bam_dt, skip_cigar = T) {
  # Quick-parse to only take reads without indels or clipping
  if (skip_cigar) {
    # Note updated CIGAR format may use X to denote mismatched bases but this is not the case
    # in our TCGA or PGP data.
    bam_dt <- bam_dt[grepl("^\\d+M$", cigar)]

    # Set rstart and rend by parsing relative to genome
    # This depends on strand orientation simply because read start and end are reversed,
    # Therefore do not apply camdac `fix_pe_strand_with_flags` before this function.
    bam_dt[, `:=`(
      rstart = (start - read.start + 1),
      rend = (end - read.start + 1)
    )]

    # Set ccgg sites to 0 for now
    bam_dt[width == 4, `:=`(rstart = 0, rend = 0)]

    return(bam_dt)
  }

  ccgg_dt <- bam_dt[width == 4]
  ccgg_dt[, `:=`(rstart = 0, rend = 0)] # Add extra columns for downstream rbind

  bam_dt <- bam_dt[width != 4]
  # Get unique reads as GAlignment object, using data.table for fast unique call
  aln <- unique(
    data.table(
      seqnames = as.character(bam_dt$chrom), pos = bam_dt$read.start, cigar = as.character(bam_dt$cigar),
      strand = GenomicAlignments::strand(bam_dt$strand), names = as.character(bam_dt$qname)
    )
  )

  aln <- GAlignments(
    seqnames = aln$seqnames, pos = aln$pos,
    cigar = aln$cigar, strand = Rle(aln$strand),
    names = aln$names
  )

  # Load SNP loci as GRanges. GenomicAlignments CIGAR parser requirement.
  gr <- GRanges(seqnames = bam_dt$chrom, ranges = IRanges(bam_dt$start, bam_dt$end))

  # Run GenomicAlignments CIGAR parser `mapToAlignments`
  res <- mapToAlignments(gr, aln)

  # Format the GRanges hits and the alignment hits as a single data table
  gh <- data.table(data.frame(gr[res$xHits]))
  names(gh) <- c("chrom", "snp_start", "snp_end", "snp_width", "snp_strand")
  ah <- data.table(data.frame(aln[res$alignmentsHits]))
  names(ah) <- gsub("^", "ah_", names(ah))
  res <- cbind(
    gh, ah, data.frame(IRanges::ranges(res)),
    data.table(qname = as.character(seqnames(res)))
  )
  rm(gh, ah, gr, aln)

  # Filter loci where the nucleotide could not be detected. These loci have
  # widths that don't match the input, which only occurs where there has been an
  # indel or soft-clipping. Propagating will only lead to NAs downstream.
  res <- res[width == snp_width]
  # Subset and rename columns
  res <- res[, .(qname, ah_strand, start, end, chrom, snp_start, snp_width)]
  names(res) <- c("qname", "strand", "rstart", "rend", "chrom", "start", "width")
  # Remove any duplicate records due to multiple hits  in `aln`.
  # Failure to do so leads to bam_dt annotated multiple timesres
  res <- unique(res)

  # Merge with BAM and return. Adds new rstart/rend columns for nucleotide
  # start and end relative to read string index
  setkey(res, qname, strand, chrom, start, width)
  setkey(bam_dt, qname, strand, chrom, start, width)
  bam_dt <- merge(bam_dt, res)

  # Return CCGG sites
  bam_dt <- rbind(bam_dt, ccgg_dt)

  return(bam_dt)
}

drop_pe_fields <- function(bam_dt) {
  pe_fields <- c("flag", "groupid", "mate_status", "cigar", "rstart", "rend")
  # Drop paired-end fields from bam_dt
  return(bam_dt[, !pe_fields, with = FALSE])
}

get_dinucs_from_seq <- function(start, read.start, seq, strand, offset = 1) {
  # Returns a vector of dinucleotide strings.
  #   Offset is required as counting and indexing differ:
  #   For example, assume the refence sequence "ATCGG" and a read aligned at "CGG",
  #   A to C is 2 steps away (start-read.start),
  #   but indexing C from the string requires position 3 (start-read.start+1)
  #   Offset differs for RRBS CCGG sites (=2) so added as argument.
  dinuc_starts <- start - read.start + offset
  dinuc_ends <- dinuc_starts + 1
  # Get dinucleotide from sequence with strand appended
  dinucs <- paste0(substr(seq, dinuc_starts, dinuc_ends), strand)
  return(dinucs)
}

get_snps_from_seq <- function(POS, read.start, seq, strand, offset = 1) {
  # Returns a vector of SNP nucleotides
  snp_position_in_seq <- POS - read.start + offset
  snp_allele <- paste0(substr(seq, snp_position_in_seq, snp_position_in_seq), strand)
  return(snp_allele)
}

get_qual_dinucs <- function(start, read.start, qual, offset = 1) {
  # Returns a vector of qual scores for the corresponding dinucleotide
  dinuc_starts <- start - read.start + offset
  dinuc_ends <- dinuc_starts + 1
  qual_dinucs <- substr(qual, dinuc_starts, dinuc_ends)
  return(qual_dinucs)
}

get_qual_snps <- function(POS, read.start, qual, offset = 1) {
  # Returns a vector of qual scores for the corresponding SNP
  snp_position_in_qual <- POS - read.start + offset
  qual_snps <- substr(qual, snp_position_in_qual, snp_position_in_qual)
  return(qual_snps)
}

get_alleles_and_qual <- function(bam_dt) {
  # Set CCGG columns (for column name continuity with RRBS version)
  bam_dt[width == 4, CCGG := TRUE]
  bam_dt[width == 4, alleles.CCGG := paste0(substr(seq, rstart, rend), strand)]
  # Set dinucleotides at CG sites for methylation rate calculation
  bam_dt[width == 2, `:=`(
    alleles.dinucs = paste0(substr(seq, rstart, rend), strand),
    qual.dinucs = paste0(substr(qual, rstart, rend), strand)
  )]

  # Set nucleotides at SNP sites, including CG-SNPs
  # As POS always gives the SNP position, we need to determine how far this is from the
  # feature position (SNP/CG/CCGG, given by `start`) and adjust the read string index `rstart` accordingly.
  # Hence, snp_pos = rstart + (POS-start)
  bam_dt[!is.na(POS), snp_idx := as.integer(rstart + (POS - start))]
  bam_dt[!is.na(POS), `:=`(
    alleles.SNP = paste0(substr(seq, snp_idx, snp_idx), strand),
    qual.SNP = paste0(substr(qual, snp_idx, snp_idx), strand)
  )]
  bam_dt[, `:=`(seq = NULL, qual = NULL, snp_idx = NULL)]
  return(bam_dt)
}

filter_bam_by_quality <- function(bam_dt, min_mapq) {
  # Set boolean flags for either SNP or dinucs above base quality 20
  # Q score encoding reference: https://support.illumina.com/help/BaseSpace_OLH_009008/Content/Source/Informatics/BS/QualityScoreEncoding_swBS.htm
  # In the regular expression, ":-@" captures all Q score ASCII symbols between ":" and "@"
  # This is capturing base quality scores at q20 and above
  bam_dt <- bam_dt[mq >= min_mapq]
  hi_qual_dinucs <- grepl("[5-9A-K:-@]{2}", bam_dt$qual.dinucs, perl = TRUE)
  hi_qual_snps   <- grepl("[5-9A-K:-@]", bam_dt$qual.SNP, perl = TRUE)

  # Set dinuc data to NA at positions where only SNP passed the quality filter
  #  This keeps the SNP data for downstream BAF/LogR but excludes the site from
  #  methylation rate calculations
  bam_dt[
    width >= 2 & !hi_qual_dinucs, 
    `:=`(alleles.dinucs = NA_character_, qual.dinucs = NA_character_)
  ]

  # Filter BAM for high quality dinucleotides and SNPs
  bam_dt <- bam_dt[
    (width >= 2 & hi_qual_dinucs) | # Hi quality dinucleotide filter
      (!is.na(POS) & hi_qual_snps) # Hi quality SNP filter
  ]

  return(bam_dt)
}

# Remove reads that do not cover entire dinucleotide at one
filter_clipped_dinucleotides <- function(bam_dt) {
  bam_dt[!(width >= 2 & read.start == end)]
}

annotate_nucleotide_counts <- function(bam_dt, rrbs = FALSE) {
  # Group everything into a single vectorized assignment block.
  bam_dt[, `:=`(
    Af  = as.integer(alleles.SNP == "A+"),
    Ar  = as.integer(alleles.SNP == "A-"),
    Cf  = as.integer(alleles.SNP == "C+"),
    Cr  = as.integer(alleles.SNP == "C-"),
    Gf  = as.integer(alleles.SNP == "G+"),
    Gr  = as.integer(alleles.SNP == "G-"),
    Tf  = as.integer(alleles.SNP == "T+"),
    Tr  = as.integer(alleles.SNP == "T-"),
    
    CGr = as.integer(alleles.dinucs == "CG-"),
    CGf = as.integer(alleles.dinucs == "CG+"),
    CAr = as.integer(alleles.dinucs == "CA-"),
    TGf = as.integer(alleles.dinucs == "TG+")
  )]

  if (rrbs) {
    bam_dt[, CCGG := as.integer(CCGG == "5pCCGG")]
  } else {
    bam_dt[, CCGG := 0L] # 0L ensures it matches the integer type of the rest
  }
  return(bam_dt)
}

flatten_pileup_to_counts <- function(bam_dt) {
  allele_counts <- bam_dt[, .(
    "Af" = sum(Af, na.rm = TRUE),
    "Ar" = sum(Ar, na.rm = TRUE),
    "Cf" = sum(Cf, na.rm = TRUE),
    "Cr" = sum(Cr, na.rm = TRUE),
    "Gf" = sum(Gf, na.rm = TRUE),
    "Gr" = sum(Gr, na.rm = TRUE),
    "Tf" = sum(Tf, na.rm = TRUE),
    "Tr" = sum(Tr, na.rm = TRUE),
    "CGf" = sum(CGf, na.rm = TRUE),
    "CGr" = sum(CGr, na.rm = TRUE),
    "TGf" = sum(TGf, na.rm = TRUE),
    "CAr" = sum(CAr, na.rm = TRUE),
    # total_depth counts reads contributing to position defined by "keyby" field,
    # hence Af selection is arbitrary and any field could be used.
    "total_depth" = .N,
    "CCGG" = sum(CCGG, na.rm = TRUE),
    "mq" = as.numeric(median(mq, na.rm = TRUE))
  ),
  keyby = .(chrom, start, end, width, POS, ref, alt)
  ]
  return(allele_counts) # df_summary in CAMDAC_RRBS
}

get_snp_allele_counts <- function(pileup_summary) {
  # Set SNP, a string of ref-alt combined from loci annotation
  pileup_summary[!is.na(ref), SNP := paste0(ref, alt)]
  # Count REF alleles. Any missing loci are set to NA
  pileup_summary[, ref_counts := data.table::fcase(
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
  pileup_summary[, alt_counts := data.table::fcase(
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
  pileup_summary[, total_counts := ref_counts + alt_counts]

  # Count all reads with nucleotides expected by CAMDAC rules. This includes
  # positions where we couldn't distinguish between bisulfite conversion and SNPs.
  # This will later be subtracted from total depth to determine unexpected nucleotide count.
  pileup_summary[, all_counts := data.table::fcase(
    is.na(ref), TGf + CAr + CGf + CGr,
    # For CT/AG SNPs, expected nucleotides are not in total_counts because they
    # confound bisulfite conversion, however we add them here for all_counts.
    SNP %like% "([GA][AG])", Af + Gf + Ar + Gr,
    SNP %like% "([CT][TC])", Cr + Tr + Tf + Cf,
    !is.na(ref), total_counts # All other positions get ref/alt counts
  )]

  # Count reads that do not have expected bases at SNP positions (i.e. SNV)
  pileup_summary[, other_counts := total_depth - all_counts]

  return(pileup_summary)
}

get_naive_snp_allele_counts <- function(pileup_summary) {
  logging::logdebug("Using naive", logger="CAMDAC")
  # Set SNP, a string of ref-alt combined from loci annotation
  pileup_summary[!is.na(ref), SNP := paste0(ref, alt)]
  # Count REF alleles. Any missing loci are set to NA
  pileup_summary[, ref_counts := data.table::fcase(
    ref == "A", Af + Ar,
    ref == "C", Cf + Cr,
    ref == "G", Gf + Gr,
    ref == "T", Tf + Tr
  )]

  # Count ALT alleles. Any loci not present are set to NA
  pileup_summary[, alt_counts := data.table::fcase(
    alt == "A", Af + Ar,
    alt == "C", Cf + Cr,
    alt == "G", Gf + Gr,
    alt == "T", Tf + Tr
  )]

  # Count total reads contributing to SNP ref/alt counts
  pileup_summary[, total_counts := ref_counts + alt_counts]

  # Count all reads with nucleotides expected by CAMDAC rules. This includes
  # positions where we couldn't distinguish between bisulfite conversion and SNPs.
  # This will later be subtracted from total depth to determine unexpected nucleotide count.
  pileup_summary[, all_counts := data.table::fcase(
    is.na(ref), TGf + CAr + CGf + CGr,
    !is.na(ref), total_counts # All other positions get ref/alt counts
  )]

  # Count reads that do not have expected bases at SNP positions (i.e. SNV)
  pileup_summary[, other_counts := total_depth - all_counts]

  return(pileup_summary)
}

# TODO: Set min_cov in CAMDAC user config
get_methylation_counts <- function(pileup_summary, min_cov) {
  # Set M from reads reporting methylation at CG dinucleotides,
  #  ignoring reads at CG-destroying SNPs.
  # Note: In RRBS version, CCGGs must be matched to SNP positions with a +1 offset
  pileup_summary[
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
  #  ignoring reads at CG-destroying SNPs
  pileup_summary[
    width > 1,
    UM := data.table::fcase(
      is.na(SNP), TGf + CAr,
      start == POS & SNP %like% "([CT][TC])", CAr,
      end == POS & SNP %like% "([AG][GA])", TGf,
      !is.na(SNP), TGf + CAr
    )
  ]

  # Set methylation counts to NA at loci with insufficient dinucleotide depth
  # This ensures that CpGs at SNP loci, where only the SNP has informative reads,
  # are not counted in downstream methylation analysis.
  pileup_summary[width > 1 & (M + UM <= min_cov), c("M", "UM") := NA]

  # Set total_counts_m as the number of reads reporting the CG methylation state
  pileup_summary[width > 1, total_counts_m := M + UM]

  return(pileup_summary)
}

get_naive_methylation_counts <- function(pileup_summary, min_cov) {
  # Set M from reads reporting methylation at CG dinucleotides,
  #  ignoring reads at CG-destroying SNPs.
  # Note: In RRBS version, CCGGs must be matched to SNP positions with a +1 offset
  pileup_summary[
    width > 1, # Calculate for CpG loci only, setting non-CpGs to NA
    M := CGf + CGr
  ]

  # Set UM as reads reporting unmethylated CG dinucleotides as above,
  #  ignoring reads at CG-destroying SNPs
  pileup_summary[
    width > 1,
    UM := TGf + CAr
  ]

  # Set methylation counts to NA at loci with insufficient dinucleotide depth
  # This ensures that CpGs at SNP loci, where only the SNP has informative reads,
  # are not counted in downstream methylation analysis.
  pileup_summary[width > 1 & (M + UM <= min_cov), c("M", "UM") := NA]

  # Set total_counts_m as the number of reads reporting the CG methylation state
  pileup_summary[width > 1, total_counts_m := M + UM]

  # Get without dinucleotide information too ----

  # Set M from reads reporting methylation at CG dinucleotides,
  #  ignoring reads at CG-destroying SNPs.
  # Note: In RRBS version, CCGGs must be matched to SNP positions with a +1 offset
  pileup_summary[
    width > 1, # Calculate for CpG loci only, setting non-CpGs to NA
    M_snuc := Cf + Cr
  ]

  # Set UM as reads reporting unmethylated CG dinucleotides as above,
  #  ignoring reads at CG-destroying SNPs
  pileup_summary[
    width > 1,
    UM_snuc := Tf + Tr
  ]

  # Set methylation counts to NA at loci with insufficient dinucleotide depth
  # This ensures that CpGs at SNP loci, where only the SNP has informative reads,
  # are not counted in downstream methylation analysis.
  pileup_summary[width > 1 & (M_snuc + UM_snuc <= min_cov), c("M_snuc", "UM_snuc") := NA]

  # Set total_counts_m as the number of reads reporting the CG methylation state
  pileup_summary[width > 1, total_counts_m_snuc := M_snuc + UM_snuc]


  return(pileup_summary)
}

filter_bad_allele_count_rows <- function(pileup_summary, min_cov) {
  # Remove positions without values for total_counts_m or total_counts
  # total_counts is NA when not a SNP. total_counts_m is NA when not met min_meth_loci read filter in previous functions.
  pileup_summary <- pileup_summary[!is.na(total_counts_m) | !is.na(total_counts)]
  loci_all <- nrow(pileup_summary) # Save initial record count for alerting users (see below)

  # Remove positions with unexpected bases at SNPs (i.e. SNV)
  pileup_summary <- pileup_summary[other_counts <= 0.05 * total_depth | (other_counts <= 1)]
  loci_low_unexpected <- nrow(pileup_summary)
  loci_high_unexpected <- loci_all - loci_low_unexpected
  # Alert users to how many loci filtered for unexpected reads
  if (loci_high_unexpected > 0) {
    # TODO: Report CpG sites filtered for logfile
  }

  # Remove positions without distinguishable ref/alt alleles or without methylation
  # Note: minimum reads filter could be implemented elsewhere in pipeline
  pileup_summary <- pileup_summary[total_counts >= min_cov | total_counts_m >= min_cov]

  return(pileup_summary)
}

compute_methylation_rates <- function(pileup_summary) {
  pileup_summary[width > 1 & total_counts_m > 0, m := M / (M + UM)]
  return(pileup_summary)
}

compute_BAFs <- function(pileup_summary) {
  pileup_summary[!is.na(ref) & total_counts > 0, BAF := alt_counts / total_counts]

  # Note: After computing BAFs, CAMDAC RRBS will remove duplicated records due
  #    to SNPs overlapping CG and CCGG sites. This is not implemented for WGBS.
  return(pileup_summary)
}

format_get_reads_result <- function(dt) {
  if (nrow(dt) == 0) {
    return(empty_count_alleles_result())
  }

  # Set a CHR column
  dt <- dt[, CHR := chrom]
  # This should always contain the CHR prefix and no sites should be NA.
  stopifnot(startsWith(as.character(dt$CHR[[1]]), "chr"))
  stopifnot(all(!is.na(dt$CHR)))

  # Format 'chrom' in Ensembl integer format for compatibility with ASCAT
  dt[, chrom := factor(
    sub("chr", "", chrom),
    levels = c(1:22, "X", "Y"),
    ordered = TRUE
  )]

  # Ensure rownames reflect row order
  rownames(dt) <- 1:nrow(dt)

  # Reorder columns
  dt <- dt[, c(
    "CHR", "chrom", "start", "end", "width", "POS", "ref", "alt", "alt_counts", "ref_counts", "total_counts", "BAF", "total_depth", "other_counts",
    "all_counts", "M", "UM", "total_counts_m", "m", "Af", "Ar", "Cf", "Cr", "Tf", "Tr", "Gf", "Gr", "CAr", "TGf", "CGf", "CGr", "CCGG", "mq"
  )]

  return(dt)
}

format_naive_get_reads_result <- function(dt) {
  if (nrow(dt) == 0) {
    return(empty_count_alleles_result())
  }

  # Set a CHR column
  dt <- dt[, CHR := chrom]
  # This should always contain the CHR prefix and no sites should be NA.
  stopifnot(startsWith(as.character(dt$CHR[[1]]), "chr"))
  stopifnot(all(!is.na(dt$CHR)))

  # Format 'chrom' in Ensembl integer format for compatibility with ASCAT
  dt[, chrom := factor(
    sub("chr", "", chrom),
    levels = c(1:22, "X", "Y"),
    ordered = TRUE
  )]

  # Ensure rownames reflect row order
  rownames(dt) <- 1:nrow(dt)

  # Reorder columns
  dt <- dt[, c(
    "CHR", "chrom", "start", "end", "width", "POS", "ref", "alt", "alt_counts", "ref_counts", "total_counts", "BAF", "total_depth", "other_counts",
    "all_counts", "M", "UM", "M_snuc", "UM_snuc", "total_counts_m", "total_counts_m_snuc", "m", "Af", "Ar", "Cf", "Cr", "Tf", "Tr", "Gf", "Gr", "CAr", "TGf", "CGf", "CGr", "CCGG", "mq"
  )]

  return(dt)
}

format_and_write_output <- function(data, output_file) {
  fs::dir_create(fs::path_dir(output_file)) # Creates directory only if it doesn't exist
  data$chrom <- as.character(data$chrom) # Ensure 'chrom' field is character - May turn to integer if X/Y regions not present
  data <- sort_genomic_dt(data) # Ensure data is sorted by chrom and POS
  data.table::fwrite(data, output_file, compress = "gzip")
  return(output_file)
}

empty_count_alleles_result <- function() {
  cols <- c(
    "CHR", "chrom", "start", "end", "width", "POS", "ref", "alt", "alt_counts", "ref_counts", "total_counts", "BAF", "total_depth", "other_counts",
    "all_counts", "M", "UM", "total_counts_m", "m", "Af", "Ar", "Cf", "Cr", "Tf", "Tr", "Gf", "Gr", "CAr", "TGf", "CGf", "CGr", "CCGG", "mq"
  )
  bam_dt <- as.data.table(matrix(NA, nrow = 0, ncol = length(cols)))
  names(bam_dt) <- cols
  return(bam_dt)
}

filter_multi_snp_loci <- function(pileup_summary) {
  # Reference files may carry potential SNPs on both CG nucleotides
  # CAMDAC currently only handles one CG-SNP pair, therefore filter for the more informative position
  # Select loci based on:
  # - One pair has a SNP, we take this one as it informs us to what degree the CG is imbalanced
  # - Take the member of the pair with the highest CpG coverage
  # - Otherwise, take the cytosine site
  # FUTURE: A tool like bis-SNP to determine SNPs per sample
  
  # Label multi_snp_loci (msl) based on duplicates
  pileup_summary[, msl := duplicated(pileup_summary, by = c("chrom", "start", "end"), fromLast = F) |
    duplicated(pileup_summary, by = c("chrom", "start", "end"), fromLast = T)]

  # Return NULL if no MSL loci found
  if(nrow(pileup_summary[msl==T])==0) {
    return(pileup_summary)
  }

  # Label loci as SNP based on BAF
  pileup_summary[msl == T, is_snp := dplyr::between(BAF, 0.1, 0.9)]

  # Label loci with maximum coverage (could be both)
  pileup_summary[msl == T, max_cov := total_counts_m == max(total_counts_m), by = c("chrom", "start", "end")]
  pileup_summary[msl == T, max_cov := ifelse(is.na(max_cov), FALSE, max_cov)]

  # Any remaining duplicates are due to equivalance in the criteria. We can therefore take the first (C) element:
  is_first <- !duplicated(pileup_summary[msl == T], by = c("chrom", "start", "end"))
  pileup_summary[msl == T, is_first := is_first]

  # Apply filters by position to determine which to keep
  pileup_summary[msl == T,
    keep := (
      (!all(is_snp) & is_snp) | # One is a SNP
        (!all(max_cov) & max_cov) # Max coverage
    ),
    by = c("chrom", "start", "end")
  ]
  pileup_summary[msl == T, keep_t := (
    (!sum(keep) == 1 & is_first) |
      (sum(keep) == 1) & keep),
  by = c("chrom", "start", "end")
  ]

  # Filter out duplicates
  pileup_summary <- pileup_summary[msl == F | (msl == T & !is.na(keep_t) & keep_t == T)]
  # Remove any remaining duplicates (potential >2 multi-SNP loci. Not currently in references.)
  pileup_summary <- pileup_summary[!duplicated(pileup_summary, by = c("chrom", "start", "end"))]

  # Remove extra columns
  pileup_summary[, `:=`(
    msl = NULL, is_snp = NULL, is_first = NULL, max_cov = NULL, keep = NULL, keep_t = NULL
  )]

  return(pileup_summary)
}

clean_and_limit_segments <- function(seg, bam_file) {
  # Set segment chromosome names based on sequence style in BAM
  # Required to read BAM correctly:
  #  e.g.: BAM having hg38 vs GRCh38 contigs 
  #  e.g.: Segments having chromosomes not currently in BAM
  bchroms <- seqnames(seqinfo(BamFile(bam_file)))
  is_ucsc <- startsWith(bchroms[[1]], "chr")
  bamstyle <- ifelse(is_ucsc, "UCSC", "NCBI")
  GenomeInfoDb::seqlevelsStyle(seg) <- bamstyle
  GenomeInfoDb::seqlevels(seg, pruning.mode="coarse") <- bchroms 
  return(seg)
}

# Wrapper ----

#' @export
cwrap_get_allele_counts <- function(bam_file, seg, loci_dt = NA, paired_end, min_mapq = 1, min_cov = 3) {
  # Loci may be NA if loci for segment chromosome (i.e. chromY)
  # is missing. Return early in these cases as no alleles to count.
  # Two conditions required to avoid error raised using is.na alone.
  # Handle the NA case from prepare_loci_list
  if (is.logical(loci_dt) && is.na(loci_dt)) {
     return(empty_count_alleles_result())
  }
  
  # Ensure loci_dt is a data.table (if it survived the NA check)
  if (!data.table::is.data.table(loci_dt)) {
     return(empty_count_alleles_result())
  }

  # Limit segments to chromosomes covered in BAM file
  seg <- clean_and_limit_segments(seg, bam_file)
  if(length(seg) == 0) {
    return(empty_count_alleles_result())
  }

  # Read BAM and annotate SNP and CPG loci
  bam_dt <- get_reads_in_segments(bam_file, seg, min_mapq, paired_end = paired_end)
  if (nrow(bam_dt) == 0) {
    return(empty_count_alleles_result())
  }

  # Overlap with loci
  bam_dt <- format_bam_for_loci_overlap(bam_dt, paired_end = paired_end)
  bam_dt <- annotate_bam_with_loci(bam_dt, loci_dt, paired_end = paired_end)
  bam_dt <- drop_positions_outside_segments(bam_dt, seg)
  if (nrow(bam_dt) == 0) {
    return(empty_count_alleles_result())
  }

  if (paired_end) {
    # For paired end samples, we must select a single read at overlapping and then fix the
    # strand column to reflect read orientation (as per single-end CAMDAC)
    bam_dt <- fix_pe_overlap_at_loci(bam_dt)
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
  bam_dt <- filter_bam_by_quality(bam_dt, min_mapq = min_mapq)
  if (nrow(bam_dt) == 0) {
    return(empty_count_alleles_result())
  }

  # Get nucleotide counts and flatten pileup
  bam_dt <- annotate_nucleotide_counts(bam_dt)
  pileup_summary <- flatten_pileup_to_counts(bam_dt)
  rm(bam_dt)

  # Apply CADMAC rules to get allele counts, methylation rates and BAFs
  pileup_summary <- get_snp_allele_counts(pileup_summary)
  pileup_summary <- get_methylation_counts(pileup_summary, min_cov)
  pileup_summary <- filter_bad_allele_count_rows(pileup_summary, min_cov)
  pileup_summary <- compute_methylation_rates(pileup_summary)
  pileup_summary <- compute_BAFs(pileup_summary)
  pileup_summary <- filter_multi_snp_loci(pileup_summary)

  # Format result for output
  result <- format_get_reads_result(pileup_summary)

  return(result)
}

