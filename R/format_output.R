##  Format allele counts output
##  Latest update: 04/05/2021
##  Version 1.0.0
##  Author: Elizabeth Larose Cadieux

#' Format output nucleotide counts
#' \code{format_output}
#'
#' @param patient_id Character variable containting the patient id number
#' @param sample_id Character variable with the sample ID
#' @param sex Character variable with  the patient expressed as "XX" for female or "XY" for male.
#' @param is_normal Logical flag set to false if the sample to be formatted is normal or tumour
#' @param path Character path variable pointing to the desired working directory.
#' This is where the output will be stored and should be constant for all CAMDAC functions.
#' Do not alter the output directory structure while running CAMDAC.
#' @param path_to_CAMDAC Character variable containting the path to the CAMDAC directory
#' including dir name  (e.g. "/path/to/CAMDAC/").
#' @param build Character variable corresponding to the reference genome used for alignment.
#' CAMDAC is compatible with "hg19", "hg38", "GRCH37","GRCH38".
#' is desired in addition to GRanges object in .RData file
#'
#' @return Concatenated SNP and CpG information

format_output <- function(patient_id, sample_id, sex,
                          is_normal = FALSE,
                          path, path_to_CAMDAC,
                          build) {
  if (getOption("scipen") == 0) {
    options(scipen = 999)
  }

  # Set output directoty
  # Do not change this - subsequent functions will look for files in this directory
  path_patient <- file.path(path, patient_id)
  path_output <- file.path(path, patient_id, "Allelecounts", sample_id)

  # Get the names of all allele counts sub-files
  index <- grepl(pattern = ".SNPs.CpGs.fst", list.files(path_output, full.names=T))
  files <- list.files(path_output, full.names=T)[index]
  rm(index)
  files <- files[order(as.numeric(gsub("^.*\\.([0-9]+)\\..*$", "\\1", files)))]

  # Combine get_allele_counts .fst outputs into one data.table
  dt_combined <- data.table(do.call(rbind, lapply(files, fst::read_fst,
    as.data.table = TRUE
  )))

  # Ensure there are no NAs before saving
  flag <- !is.na(dt_combined$CHR)
  if (sum(flag) < length(flag)) {
    dt_combined[, ] <- dt_combined[flag, ]
  }
  rm(flag)

  # Ensure that spurious alignments to Y in females are removed
  if (sex == "XX") {
    dt_combined <- dt_combined[!dt_combined$chrom %in% c("Y", "chrY"), ]
  }

  # Create allele counts combined output file name
  output_file_prefix <- paste0(path_output, "/", patient_id, ".", sample_id)
  f_nm <- paste0(output_file_prefix, ".SNPs.CpGs.all.sorted.RData")

  # Save tumour data
  if (is_normal == FALSE) {
    dt_tumour <- dt_combined
    save(dt_tumour, file = f_nm)
  }

  # Save normal data
  if (is_normal == TRUE) {
    dt_normal <- dt_combined

    # save
    save(dt_normal, file = f_nm)
  }

  # clean up
  if (file.exists(f_nm)) {
    file.remove(files)
  }
  rm(f_nm, files)

  # Create dir and set output file for msp1 fragments sizes and copy number
  dir.create(file.path(path_patient, "Copy_number", sample_id), recursive = TRUE)
  outfile_prefix <- paste0(path_patient, "/Copy_number/", sample_id, "/")

  # Get normal msp1 fragments sizes and nucleotide content
  get_msp1_fragments(
    dt = dt_combined, build = build,
    path_to_CAMDAC = path_to_CAMDAC,
    outfile = outfile_prefix
  )
  cat("Msp1 fragments information obtained for patient\n")
  rm(outfile_prefix)
}


#' @title get_msp1_fragments
#' @description get msp1 fragments
#' @param dt data.table object with containing all covered CCGGs in the sample
#' @param build Character, Either "hg19", "hg38", "GRCH37","GRCH38"
#' @param path_to_CAMDAC Character string containting the path to the CAMDAC dir including
#' dir name e.g. "~/CAMDAC/"
#' @param outfile character srting with output filename
#'
#' @author elizabeth larose cadieux

get_msp1_fragments <- function(dt, build, path_to_CAMDAC, outfile) {
  # Set build to to assembly version disregarging USCS vs. Ensembl
  if (build == "GRCH37") {
    build <- "hg19"
  }
  if (build == "GRCH38") {
    build <- "hg38"
  }

  # Load fragments file
  msp1_fragments_file <- file.path(
    path_to_CAMDAC,
    paste0("pipeline_files/msp1_fragments/msp1_fragments_RRBS_", build, ".fst")
  )
  fragments <- fst::read_fst(path = msp1_fragments_file, as.data.table = TRUE)

  # Assign CpG IDs
  dt = data.table::data.table(dt)
  dt[, CpG_ID := paste(CHR, start, end, sep = "_")]
  dt <- dt[!duplicated(CpG_ID), ]

  # Format reference fragments file to correct annotation set
  ch <- dt[1, as.character(CHR)]
  levs <- paste0("chr", c(1:22, "X", "Y"))
  if (grepl("chr", ch) == FALSE) {
    dt[, CHR := paste0("chr", CHR)]
  }
  rm(ch)

  # Get CCGGs with coverage in the sample
  dt_CCGGs <- dt[total_counts_m > 0 & CCGG > 0, ]
  dt_CCGGs <- dt_CCGGs[, .SD, .SDcols = c("CHR", "start", "end")]
  setkeyv(dt_CCGGs, cols = c("CHR", "start", "end"))

  # Get all positions with coverage in the sample
  dt <- dt[, .SD, .SDcols = c("CHR", "start", "end", "CpG_ID")]
  setkeyv(dt, cols = c("CHR", "start", "end"))

  # Assign fragments IDs
  fragments[, ID := 1:nrow(fragments)]

  # Format reference fragments file to correct annotation set
  ch <- fragments[1, as.character(CHR)]
  if (grepl("chr", ch) == FALSE) {
    fragments[, CHR := paste0("chr", CHR)]
  }
  rm(ch)
  setkeyv(fragments, cols = c("CHR", "start", "end"))

  # Each reference 5'CCGG (including common polymorphisms)
  # should have read density if a CCGG is present in the data
  fragments_patient <- foverlaps(dt_CCGGs, fragments, nomatch = 0)
  rm(dt_CCGGs)
  fragments_patient[, c("i.start", "i.end") := NULL]
  setkeyv(fragments_patient, cols = c("CHR", "start", "end"))

  # Get CpGs on mspI fragments of known length
  overlaps <- foverlaps(dt, fragments_patient, nomatch = 0)
  overlaps[, c("i.start", "i.end") := NULL]

  # In some cases, only one of the two mspI fargment 5'CCGG will have density
  # (either by chance at low coverage sites, due to one of the ends of the
  # fragment having low mappability or, less likely, due to SVs/SNVs forming new fragment ends).
  # This problem is observed ~ 1% of RRBS loci.

  # Get CpGs on fragments of unknown length
  unknown_fragments <- dt[!CpG_ID %in% unique(overlaps$CpG_ID), ]
  # nrow(unknown_fragments) / nrow(dt)
  setkeyv(unknown_fragments, cols = c("CHR", "start", "end"))

  # Add the in sillico fragments for those where only one 5'end has reads
  overlaps <- foverlaps(unknown_fragments, fragments, nomatch = 0)
  rm(fragments)
  overlaps[, c("i.start", "i.end", "CpG_ID") := NULL]

  # Remove loci for which we do not have a good approx in silico fragment estimate
  overlaps <- overlaps[msp1_length <= 1000, ]

  # Concatenate both sets
  fragments_patient <- rbind(fragments_patient, overlaps)
  rm(overlaps)
  fragments_patient <- fragments_patient[!duplicated(fragments_patient)]

  # Order data.table
  fragments_patient <-
    fragments_patient[order(factor(CHR, levels = levs, ordered = TRUE), start, end)]

  # Remove loci for which we do not have a good approx in silico fragment estimate
  # fragments_patient[, table(msp1_length<=1000)] # This steps removes ~ 5% of CpG loci
  fragments_patient <- fragments_patient[msp1_length <= 1000, ]

  # First position of each fragment will always be cleaved by MspI digestion
  # Leaving a 5' CGG on the (+) strand and CCG on the (-) strand in single end directional RRBS
  fragments_patient[, c("start", "end") := .(start + 1, end - 1)]

  # Save fragments file
  save(fragments_patient, file = paste0(outfile, "msp1_fragments_RRBS.RData"))

  # Extract fragment size and log10 transformed values
  df_fragments <- fragments_patient[, .(l = msp1_length, ll = log10(msp1_length))]

  # Remove fragments smaller than 35bp before calculating mean statistics
  df_fragments <- df_fragments[l > 35, ]

  # Get mean fragment length and inter-quartile range
  mean_length <- df_fragments[, round2(mean(l, na.rm = TRUE), digits = 0)]
  q25 <- df_fragments[, quantile(l, 0.25)]
  q75 <- df_fragments[, quantile(l, 0.75)]

  # plot log10 fragment size distribution
  outfile <- paste0(outfile, "fragment_length_histogram.pdf")
  p <- ggplot2::ggplot(df_fragments) +
    geom_histogram(aes(x = l, y = ..count..), col = "cornflowerblue", fill = "white", bins = 100) +
    theme_classic() +
    ylab("Number of fragments") +
    coord_cartesian(xlim = c(35, 1000)) + #+coord_cartesian(xlim=c(log10(40),3))+
    # scale_x_continuous(name="Log Msp1 fragment length", breaks = seq(1,3,.25),labels = round2(10^seq(1,3,.25), digits = 0))+
    scale_x_continuous(name = "MspI fragment length", breaks = seq(100, 1000, 100), labels = seq(100, 1000, 100)) +
    ggtitle(paste0(
      "Mean MspI fragment length = ", mean_length,
      "bp\nand the inter-quartile range is [", q25, ", ", q75, "] bp"
    )) +
    theme(plot.title = element_text(hjust = 0.5))
  ggsave(file = outfile, plot = p, device = "pdf", width = 4.5, height = 3, units = "in")
}


#' @description Round numerical values to 'n' digits
#' @param x Numerical vector containing the numbers to round
#' @param digits Numerical value representing the number of decimal digits to retain
#' @return rounded numerical vector
round2 <- function(x, digits) {
  ifelse(as.integer(x * (10^(digits + 1))) %% 10 >= 5, ceiling(x * (10^digits)) / (10^digits), floor(x * (10^digits)) / (10^digits))
}

# END
