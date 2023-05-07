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

    # Select essential fields
    ab = ol[, .(major_cn, minor_cn, hap_BAF, hap_reads)]
    ab$cnstate = paste0(ab$major_cn, "+", ab$minor_cn)
    ab$cnix = seq_len(nrow(ab))
    
    # Split by balance
    ab_bal = ab[ major_cn == minor_cn,  ]
    ab_imbal = ab[ major_cn != minor_cn,  ]

    # Assign major to ref or alt allele
    ab_bal[, maj_assign:=maj_by_gauss(hap_BAF, balanced=T)$maj, by=cnstate]
    ab_imbal[, maj_assign:=maj_by_gauss(hap_BAF, balanced=F)$maj, by=cnstate]
    ab_phas = rbind(ab_bal, ab_imbal, ab[is.na(major_cn) | is.na(minor_cn)], fill=T)[order(cnix), ]

    # Assign ref and alt CN
    ab_phas[, ref_CN:=data.table::fcase(
        maj_assign == "ref", major_cn,
        maj_assign == "alt", minor_cn,
        maj_assign == "balanced", major_cn
    )]
    ab_phas[, alt_CN:=data.table::fcase(
        maj_assign == "ref", minor_cn,
        maj_assign == "alt", major_cn,
        maj_assign == "balanced", major_cn
    )]

    # Cleanup fields
    ab_phas$cnix=NULL
    add_cols = setdiff(names(ab_phas), names(ol))
    res = cbind(ol, ab_phas[, ..add_cols])
    return(res)
}

maj_by_gauss <- function(BAF, balanced=T){
    # Reserve original BAF values to determine ref or alt as major or minor
    BAF_orig = BAF
    # First, ensure BAF is mirrored by flipping to 0.5 and 1 window
    BAF[ BAF < 0.5 ] = 1 - BAF[ BAF < 0.5 ]

    # Classify BAF as major, minor, balanced, or NA
      # truncate data to values less than 1. 
    truncated_data <- BAF[BAF != 1]
    if(length(truncated_data) == 0){
        return(list("maj" = "NA", "params" = c("NA", "NA", "NA", "NA")))
    }

    # Get thresholds from gaussian fit. Does not fit to extreme values on 0 or 1
    # This avoids fitting to outliers yet should identify peak regardless
    fit <- fit_gaussian(truncated_data)
    thresh = qnorm(c(0.025, 0.975), fit$mean, fit$sd)
    lower = thresh[[1]]
    upper = thresh[[2]]
    
    # Set default threshold for filtering out hets in imbalances
    hets_upper = qnorm(c(0.99), mean=0.5, sd=0.01)

    # If state is balanced, return "balanced"
    if(balanced){
        maj = ifelse(BAF < upper, "balanced", "NA")
    }else{
        # One model fits at peak, another fits at 0.5
        in_dist = dplyr::between(BAF, lower, upper)
        is_ref = BAF_orig < 0.5
        is_het_outdist = BAF < hets_upper
        maj = dplyr::case_when(
            is_het_outdist ~ "NA",
            is_ref & in_dist ~ "ref",
            !is_ref & in_dist ~ "alt",
            TRUE ~ "NA"
        )
    }

    # Setup parameters for return
    params = c(
        mean = fit$mean,
        sd = fit$sd,
        qlower = lower,
        qupper = upper
    )
    
    return(list("maj" = maj, "params" = params))
}

fit_gaussian <- function(data) {
  fit <- fitdistr(data, "normal") # fit Gaussian distribution to truncated data
  mean <- fit$estimate[1] # extract mean from fitted parameters
  sd <- fit$estimate[2] # extract standard deviation from fitted parameters
  return(list(mean = mean, sd = sd))
}
