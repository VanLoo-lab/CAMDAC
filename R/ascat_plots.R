
#' @title ascat.m.plotRawData
#' @description Plot tumour and germline BAF and LogR
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#' @param raw_LogR vector with the LogR values before correction # NM: Not used so removed
#' @param pch type of data points in plot
#' @param cex size of data points in plot
#' @param lim_logR y-axis limits on logR plot
#'
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#' @author Peter Van Loo
#' @noRd
ascat.m.plotRawData <- function(ASCATobj, outdir, pch = 10, cex = 0.2, lim_logR = 2.5) {
  print.noquote("Plotting tumor data")
  for (i in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    colls <- ifelse(ASCATobj$Germline_BAF[, i] < 0.85 & ASCATobj$Germline_BAF[, i] > 0.15, "red", "red")
    # set point colours to show SNP germline genotype
    outfile_t <- fs::path(outdir, paste(ASCATobj$samples[i], ".tumour.png", sep = ""))
    png(filename = outfile_t, width = 2000, height = 1250, res = 200)
    par(
      mar = c(0.5, 5, 5, 0.5), mfrow = c(3, 1), cex = 0.4, cex.main = 3, cex.axis = 2,
      pch = ifelse(dim(ASCATobj$Tumor_LogR)[1] > 100000, ".", 20)
    )
    plot(c(1, dim(ASCATobj$Tumor_LogR)[1]), c(-lim_logR, lim_logR),
      type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, LogR", sep = ""),
      xlab = "", ylab = ""
    )
    points(ASCATobj$Tumor_LogR[, i], col = "red", cex = 0.2)
    # points(ASCATobj$Tumor_LogR[,i],col=rainbow(24)[ASCATobj$SNPpos$Chr])
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk <- ASCATobj$ch[[j]]
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 2, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    # Note: no corrected data currently passed NM edit
    # plot(c(1,dim(ASCATobj$Tumor_LogR)[1]), c(-lim_logR ,lim_logR ),
    #      type = "n", xaxt = "n", main = paste(ASCATobj$samples[i], ", tumor data, corrected LogR", sep = ""),
    #      xlab = "", ylab = "")
    # points(ASCATobj$Tumor_LogR[,i],col="red")
    # #points(ASCATobj$Tumor_LogR[,i],col=rainbow(24)[ASCATobj$SNPpos$Chr])
    # abline(v=0.5,lty=1,col="lightgrey")
    # chrk_tot_len = 0
    # for (j in 1:length(ASCATobj$ch)) {
    #   chrk = ASCATobj$ch[[j]];
    #   chrk_tot_len_prev = chrk_tot_len
    #   chrk_tot_len = chrk_tot_len + length(chrk)
    #   vpos = chrk_tot_len;
    #   tpos = (chrk_tot_len+chrk_tot_len_prev)/2;
    #   text(tpos,2,ASCATobj$chrs[j], pos = 1, cex = 2)
    #   abline(v=vpos+0.5,lty=1,col="lightgrey")
    # }
    plot(c(1, dim(ASCATobj$Tumor_BAF)[1]), c(0, 1),
      type = "n", xaxt = "n",
      main = paste(ASCATobj$samples[i], ", tumor data, BAF", sep = ""), xlab = "", ylab = ""
    )
    points(ASCATobj$Tumor_BAF[, i], col = colls, pch = pch, cex = 0.2)
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk <- ASCATobj$ch[[j]]
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    dev.off()
  }

  if (!is.null(ASCATobj$Germline_LogR)) {
    print.noquote("Plotting germline data")
    for (i in 1:dim(ASCATobj$Germline_LogR)[2]) {
      outfile_g <- fs::path(outdir, paste(ASCATobj$samples[i], ".germline.png", sep = ""))
      png(filename = outfile_g, width = 2000, height = 750, res = 200)
      par(
        mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4, cex.main = 3, cex.axis = 2,
        pch = ifelse(dim(ASCATobj$Tumor_LogR)[1] > 100000, ".", 20)
      )
      plot(c(1, dim(ASCATobj$Germline_LogR)[1]), c(-1, 1),
        type = "n", xaxt = "n",
        main = paste(ASCATobj$samples[i], ", germline data, LogR", sep = ""), xlab = "", ylab = ""
      )
      points(ASCATobj$Germline_LogR[, i], col = "red", cex = 0.2)
      abline(v = 0.5, lty = 1, col = "lightgrey")
      chrk_tot_len <- 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk <- ASCATobj$ch[[j]]
        chrk_tot_len_prev <- chrk_tot_len
        chrk_tot_len <- chrk_tot_len + length(chrk)
        vpos <- chrk_tot_len
        tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
        text(tpos, 2, ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
      }
      plot(c(1, dim(ASCATobj$Germline_BAF)[1]), c(0, 1),
        type = "n", xaxt = "n",
        main = paste(ASCATobj$samples[i], ", germline data, BAF", sep = ""), xlab = "", ylab = ""
      )
      points(ASCATobj$Germline_BAF[, i], col = colls, pch = pch, cex = 0.2)
      abline(v = 0.5, lty = 1, col = "lightgrey")
      chrk_tot_len <- 0
      for (j in 1:length(ASCATobj$ch)) {
        chrk <- ASCATobj$ch[[j]]
        chrk_tot_len_prev <- chrk_tot_len
        chrk_tot_len <- chrk_tot_len + length(chrk)
        vpos <- chrk_tot_len
        tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
        text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
        abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
      }
      dev.off()
    }
  }
}


#' @title ascat.m.plotSegmentedData
#' @description Plot segmentated BAF LogR
#' @param ASCATobj an ASCAT object (e.g. data structure from ascat.loadData)
#'
#' @return Produces png files showing the logR and BAF values for tumour and germline samples
#' @author Peter Van Loo
#' @noRd
ascat.m.plotSegmentedData <- function(ASCATobj, fname = "", outdir, lim_logR = 2) {
  for (arraynr in 1:dim(ASCATobj$Tumor_LogR)[2]) {
    outfile_t <- fs::path(outdir, paste(ASCATobj$samples[arraynr], ".ASPCF.png",
      sep = ""
    ))
    Select_nonNAs <- rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]])
    AllIDs <- 1:dim(ASCATobj$Tumor_LogR)[1]
    names(AllIDs) <- rownames(ASCATobj$Tumor_LogR)
    HetIDs <- AllIDs[Select_nonNAs]
    png(filename = outfile_t, width = 2000, height = 1000, res = 200)
    par(
      mar = c(0.5, 5, 5, 0.5), mfrow = c(2, 1), cex = 0.4,
      cex.main = 3, cex.axis = 2
    )
    r <- ASCATobj$Tumor_LogR_segmented[
      rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]),
      arraynr
    ]
    beta <- ASCATobj$Tumor_BAF_segmented[[arraynr]][, , drop = FALSE]
    plot(c(1, length(r)), c(-lim_logR, lim_logR),
      type = "n", xaxt = "n",
      main = paste(fname, ", LogR", sep = ""), xlab = "", ylab = ""
    )
    points(ASCATobj$Tumor_LogR[
      rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]),
      arraynr
    ], col = rgb(1, 0, 0, 0.5), pch = ".", cex = 0.20)
    points(r, col = "blue", cex = 0.2)
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk <- intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, lim_logR - 0.5, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    plot(c(1, length(beta)), c(0, 1),
      type = "n", xaxt = "n",
      main = paste(fname, ", BAF", sep = ""), xlab = "", ylab = ""
    )
    points(ASCATobj$Tumor_BAF[
      rownames(ASCATobj$Tumor_BAF_segmented[[arraynr]]),
      arraynr
    ], col = rgb(1, 0, 0, 0.5), pch = ".", cex = 0.20)
    points(beta, col = "blue", cex = 0.2)
    points(1 - beta, col = "blue", cex = 0.2)
    abline(v = 0.5, lty = 1, col = "lightgrey")
    chrk_tot_len <- 0
    for (j in 1:length(ASCATobj$ch)) {
      chrk <- intersect(ASCATobj$ch[[j]], HetIDs)
      chrk_tot_len_prev <- chrk_tot_len
      chrk_tot_len <- chrk_tot_len + length(chrk)
      vpos <- chrk_tot_len
      tpos <- (chrk_tot_len + chrk_tot_len_prev) / 2
      text(tpos, 1, ASCATobj$chrs[j], pos = 1, cex = 2)
      abline(v = vpos + 0.5, lty = 1, col = "lightgrey")
    }
    dev.off()
  }
}

# Copied and adjusted from Tom's github repo:
# https://github.com/tlesluyes/ascat/blob/dc53b739504be83d309be6cffa62c16a1de770df/ASCAT/R/ascat.plots.R#L389
ascat.m.plotAdjustedAscatProfile <- function(camdac_cna, outfile_name, sample_id = "SAMPLE", REF = "hg38", y_limit = 5, plot_unrounded = F, battenberg = F, ascat_colours = F) {
  # convert CAMDAC_CNA to expected values
  ASCAT_output_object <- list()
  ASCAT_output_object$segments_raw <- data.frame(
    camdac_cna$cna[, .(sample = sample_id, chr = chrom, startpos = start, endpos = end, nMajor = nA, nMinor = nB, nAraw = nA, nBraw = nB)]
  )
  ASCAT_output_object$segments <- data.frame(
    camdac_cna$cna[, .(sample = sample_id, chr = chrom, startpos = start, endpos = end, nMajor = nA, nMinor = nB)]
  )
  # Set variables for final section of plot
  SAMPLE <- sample_id
  ASCAT_output_object$purity[SAMPLE] <- camdac_cna$purity
  ASCAT_output_object$ploidy[SAMPLE] <- camdac_cna$ploidy
  ASCAT_output_object$goodnessOfFit[SAMPLE] <- camdac_cna$fit
  ASCAT_output_object$nonaberrantarrays[SAMPLE] <- F

  if (plot_unrounded) {
    SEGMENTS <- ASCAT_output_object$segments_raw[, c(1:4, 7:8)]
    colnames(SEGMENTS)[5:6] <- c("nMajor", "nMinor")
    SEGMENTS$nMajor <- SEGMENTS$nMajor + SEGMENTS$nMinor
    colourA <- "#c725e3" # purple
    colourB <- "#e37825" # orange
  } else {
    SEGMENTS <- ASCAT_output_object$segments
    SEGMENTS$nMajor <- SEGMENTS$nMajor - 0.1
    SEGMENTS$nMinor <- SEGMENTS$nMinor + 0.1
    colourA <- "#e03546" # red
    colourB <- "#3557e0" # blue
  }
  SEGMENTS$nMajor <- ifelse(SEGMENTS$nMajor > y_limit, y_limit + 0.1, SEGMENTS$nMajor)
  SEGMENTS$nMinor <- ifelse(SEGMENTS$nMinor > y_limit, y_limit + 0.1, SEGMENTS$nMinor)

  if (battenberg) {
    colourA <- "#e4a329"
    colourB <- "#000000"
  }
  if (ascat_colours) {
    colourA <- "#00fd31"
    colourB <- "#fd2b1a"
  }

  if (REF == "hg19") {
    REF <- data.frame(
      chrom = c(1:22, "X"),
      start = rep(1, 23),
      end = c(
        249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431,
        135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248,
        59128983, 63025520, 48129895, 51304566, 155270560
      )
    )
  } else if (REF == "hg38") {
    REF <- data.frame(
      chrom = c(1:22, "X"),
      start = rep(1, 23),
      end = c(
        248956422, 242193529, 198295559, 190214555, 181538259, 170805979, 159345973, 145138636, 138394717,
        133797422, 135086622, 133275309, 114364328, 107043718, 101991189, 90338345, 83257441, 80373285,
        58617616, 64444167, 46709983, 50818468, 156040895
      )
    )
  } else {
    stopifnot(is.data.frame(REF))
    stopifnot(identical(colnames(REF), c("chrom", "start", "end")))
  }

  SEGMENTS$chr <- gsub("^chr", "", SEGMENTS$chr)
  stopifnot(all(ASCAT_output_object$segments$chr %in% REF$chrom))
  REF$size <- REF$end - REF$start + 1
  REF$middle <- 0
  for (i in 1:nrow(REF)) {
    if (i == 1) {
      REF$middle[i] <- REF$size[i] / 2
    } else {
      REF$middle[i] <- sum(as.numeric(REF$size[1:(i - 1)])) + REF$size[i] / 2
    }
  }
  rm(i)
  REF$cumul <- cumsum(as.numeric(REF$size))
  REF$add <- cumsum(as.numeric(c(0, REF$size[1:(nrow(REF) - 1)])))

  SEGMENTS$startpos_adjusted <- SEGMENTS$startpos
  SEGMENTS$endpos_adjusted <- SEGMENTS$endpos
  for (CHR in unique(REF$chrom)) {
    INDEX <- which(SEGMENTS$chr == CHR)
    if (length(INDEX) > 0) {
      SEGMENTS$startpos_adjusted[INDEX] <- SEGMENTS$startpos_adjusted[INDEX] + REF$add[which(REF$chrom == CHR)]
      SEGMENTS$endpos_adjusted[INDEX] <- SEGMENTS$endpos_adjusted[INDEX] + REF$add[which(REF$chrom == CHR)]
    }
    rm(INDEX)
  }
  rm(CHR)


  for (SAMPLE in sort(unique(SEGMENTS$sample))) {
    SEGS <- SEGMENTS[which(SEGMENTS$sample == SAMPLE), ]
    if (nrow(SEGS) == 0) warning(paste0("No segments for sample: ", SAMPLE))
    maintitle <- paste(SAMPLE, "  Ploidy: ", sprintf("%1.2f", ASCAT_output_object$ploidy[SAMPLE]), ", purity: ", sprintf("%2.0f", ASCAT_output_object$purity[SAMPLE] * 100), "%, goodness of fit: ", sprintf("%2.1f", ASCAT_output_object$goodnessOfFit[SAMPLE]), "%", ifelse(isTRUE(ASCAT_output_object$nonaberrantarrays[SAMPLE]), ", non-aberrant", ""), sep = "")
    png(filename = outfile_name, width = 2000, height = (y_limit * 100), res = 200)
    par(mar = c(0.5, 5, 5, 0.5), cex = 0.4, cex.main = 3, cex.axis = 2.5)
    ticks <- seq(0, y_limit, 1)
    plot(c(1, REF$cumul[nrow(REF)]), c(0, y_limit), type = "n", xaxt = "n", yaxt = "n", main = maintitle, xlab = "", ylab = "")
    axis(side = 2, at = ticks)
    abline(h = ticks, col = "lightgrey", lty = 1)
    rect(SEGS$startpos_adjusted, (SEGS$nMajor - 0.07), SEGS$endpos_adjusted, (SEGS$nMajor + 0.07), col = ifelse(SEGS$nMajor >= y_limit, adjustcolor(colourA, red.f = 0.75, green.f = 0.75, blue.f = 0.75), colourA), border = ifelse(SEGS$nMajor >= y_limit, adjustcolor(colourA, red.f = 0.75, green.f = 0.75, blue.f = 0.75), colourA))
    rect(SEGS$startpos_adjusted, (SEGS$nMinor - 0.07), SEGS$endpos_adjusted, (SEGS$nMinor + 0.07), col = ifelse(SEGS$nMinor >= y_limit, adjustcolor(colourB, red.f = 0.75, green.f = 0.75, blue.f = 0.75), colourB), border = ifelse(SEGS$nMinor >= y_limit, adjustcolor(colourB, red.f = 0.75, green.f = 0.75, blue.f = 0.75), colourB))
    abline(v = c(1, REF$cumul), lty = 1, col = "lightgrey")
    text(REF$middle, y_limit, REF$chrom, pos = 1, cex = 2)
    dev.off()
    rm(SEGS, ticks, maintitle)
  }
  rm(SAMPLE)
}
