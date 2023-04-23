cmain_asm_allele_counts <- function(sample, config) {
    loginfo("ASM allele counting for %s", paste0(sample$patient_id, ":", sample$id))

    # Skip if no BAM file provided
    if (is.null(sample$bam)) {
        loginfo("No BAM. Skipping allele counting for %s", paste0(sample$patient_id, ":", sample$id))
        return(NULL)
    }

    #  Check if outputs exist and skip if required
    outfile <- get_fpath(sample, config, "asm_counts")
    if (file.exists(outfile) && !config$overwrite) {
        logwarn("Skipping ASM allele counting for %s", paste0(sample$patient_id, ":", sample$id))
        return(outfile)
    }

    # Create temporary directory for allele counts files
    # Use tempfile to create unique suffix and avoid overwrites on failed runs
    tempdir <- tempfile(pattern = "asm_counts", tmpdir = get_fpath(sample, config, "asm_counts", dir = T))
    fs::dir_create(tempdir)

    # Get SNP loci as segments to analyse. Parallelised over config$n_seg_split
    snps_gr <- load_asm_snps_gr(sample, config)
    # If regions to analyse are given, limit SNPs to regions
    if (!is.null(config$regions)) {
        segments <- read_segments_bed(config$regions)
        seg_gr <- Reduce(c, segments)
        snps_gr <- subsetByOverlaps(snps_gr, seg_gr)
    }

    # Split by chromosome for parallelisation
    snps_grl <- split(snps_gr, seqnames(snps_gr))


    # List files containing SNP and CpG loci for reference genome
    loci_files <- get_reference_files(config, type = "loci_files")

    # Load sample data
    bam_file <- sample$bam
    paired_end <- is_pe(config)
    drop_ccgg <- is_ccgg(config)
    min_mapq <- config$min_mapq
    min_cov <- config$min_cov

    # Initialise parallel workers.
    doParallel::registerDoParallel(cores = config$n_cores)

    #   Set warn=2 to ensure foreach fails if any of the parallel workers are terminated due to memory.
    #   without this option, foreach simply returns a warning and software continues
    options(warn = 2)

    tmpfiles <- foreach(seg = snps_grl, .combine = "c") %dopar% {
        # Loop over SNPs to phase
        loci_dt <- load_asm_loci_for_segment(seg, loci_files)
        ac_file <- cwrap_asm_get_allele_counts(bam_file, seg, loci_dt, paired_end, drop_ccgg, min_mapq = min_mapq, min_cov = min_cov)
        tmp <- tempfile(tmpdir = tempdir, fileext = ".qs")
        qs::qsave(ac_file, tmp)
        rm(loci_dt, ac_file, seg)
        gc()
        return(tmp)
    }
    options(warn = 0)

    # Define function to combine the allele counts objects into a single list
    bind_asm_obs <- function(x, y) {
        nobj <- list()
        nobj$asm_cg <- rbind(x$asm_cg, y$asm_cg, fill = T)
        nobj$hap_stats <- rbind(x$hap_stats, y$hap_stats, fill = T)
        nobj$map <- rbind(x$map, y$map, fill = T)
        return(nobj)
    }
    # Combine temporary files with allele counts results into a single data table
    result <- foreach(i = tmpfiles, .combine = bind_asm_obs) %dopar% {
        qs::qread(i)
    }

    # Write to output(s) file
    asm_ac_out <- write_asm_counts_output(result, sample, config)

    # Delete temporary files
    fs::dir_delete(tempdir)

    # Stop parallel workers. When running the pipeline multiple times in an R session,
    # R re-uses workers but does not clear memory. Hence large objects in foreach loops will remain.
    doParallel::stopImplicitCluster()
    return(asm_ac_out)
}

cmain_asm_make_methylation <- function(sample, config) {
    # Skip if asm_counts_file doesn't exist
    asm_counts_file <- get_fpath(sample, config, "asm_counts")
    if (!file.exists(asm_counts_file)) {
        logwarn("No ASM allele counts. Skipping ASM methylation for %s", paste0(sample$patient_id, ":", sample$id))
        stop()
        return(NULL)
    }

    # Skip if asm_meth_file exists
    asm_meth_file <- get_fpath(sample, config, "asm_meth")
    if (file.exists(asm_meth_file) && !config$overwrite) {
        logwarn("ASM methylation already exists. Skipping ASM methylation for %s", paste0(sample$patient_id, ":", sample$id))
        return(asm_meth_file)
    }

    # Load DNA methylation object for asm
    asm_counts <- fread_chrom(asm_counts_file)

    # Select DNA methylation fields
    asm_meth <- asm_counts[
        width == 2,
        .(chrom, start, end, alt_total_counts_m, ref_total_counts_m, alt_m, ref_m)
    ]

    # Save ASM methylation
    asm_meth_outfile <- get_fpath(sample, config, "asm_meth")
    fs::dir_create(fs::path_dir(asm_meth_outfile))
    data.table::fwrite(asm_meth, asm_meth_outfile)

    return(asm_meth_outfile)
}


cmain_asm_make_snps <- function(tumor, germline, infiltrates, origin, config) {
    # Check that ASM snps file is availabe for tumor. If so, return NULL
    asm_snps_file <- get_fpath(tumor, config, "asm_snps")
    if (file.exists(asm_snps_file)) {
        loginfo("ASM snps file found for tumor.")
        # Loop over remaining objects
        for (i in list(germline, infiltrates, origin)) {
            if (is.null(i)) {
                next
            }
            # Add ASM snps from tumor to object if currently Null
            i_asm_snps <- get_fpath(i, config, "asm_snps")
            if (!file.exists(i_asm_snps)) {
                loginfo("Attaching existing ASM SNPs to %s", i$id)
                attach_output(i, config, "asm_snps", asm_snps_file)
            }
        }
        return(NULL)
    }

    # Else, raise an error if germline is NULL
    if (is.null(germline)) {
        stop(paste0(
            "No ASM snps file available for tumor. Germline sample is required to extract SNPs from,",
            "or ASM snps file must be provided for tumor object."
        ))
    }

    # If ASM snps are not available for the tumor, run bulk allele-counts on germline and extract SNPs
    loginfo("ASM snps file not found for tumor. Extracting SNPs from germline for tumor.")
    cmain_count_alleles(germline, config)
    cmain_make_snps(germline, config)

    # Load het SNPs
    nsnps_f <- get_fpath(germline, config, "snps")
    n_snp <- fread_chrom(nsnps_f)
    n_snp <- n_snp[dplyr::between(BAF, 0.1, 0.9), .(chrom, pos = POS, ref, alt, BAF)]

    # Save hets as ASM SNPs for germline
    n_asm_snps_file <- get_fpath(germline, config, "asm_snps")
    fs::dir_create(fs::path_dir(n_asm_snps_file))
    data.table::fwrite(n_snp, n_asm_snps_file)

    # Save hets as ASM SNPs for tumor
    fs::dir_create(fs::path_dir(asm_snps_file))
    data.table::fwrite(n_snp, asm_snps_file)

    # Save hets as ASM SNPs for origin and infiltrates if present
    if (!is.null(infiltrates)) {
        i_asm_snps_file <- get_fpath(infiltrates, config, "asm_snps")
        fs::dir_create(fs::path_dir(i_asm_snps_file))
        data.table::fwrite(n_snp, i_asm_snps_file)
    }
    if (!is.null(origin)) {
        o_asm_snps_file <- get_fpath(origin, config, "asm_snps")
        fs::dir_create(fs::path_dir(o_asm_snps_file))
        data.table::fwrite(n_snp, o_asm_snps_file)
    }
    loginfo("ASM SNPS file created from germline for: {tumor$patient_id}:{tumor$id}")
}

cmain_asm_call_cna <- function(tumor, germline, config) {
    # Check that CNA file is available for the tumor
    asm_cna_file <- get_fpath(tumor, config, "asm_cna")
    if (fs::file_exists(asm_cna_file)) {
        loginfo("CNA file found for tumor.")
        return(NULL)
    }

    # Preprocess CpG, SNP and methylation data for all samples
    preprocess(
        list(tumor, germline),
        config
    )

    # Combine tumor-germline SNPs and call CNAs
    cmain_bind_snps(tumor, germline, config)
    cmain_call_cna(tumor, germline, config)
    cna_file <- get_fpath(tumor, config, "cna")
    attach_output(tumor, config, "asm_cna", cna_file)
    loginfo("CNA file created for tumor: {tumor$patient_id}:{tumor$id}")
}

cmain_fit_meth_cna <- function(tumor, config) {
    # Skip if meth_cn file exists
    asm_meth_cna <- get_fpath(tumor, config, "asm_meth_cna")
    if (file.exists(asm_meth_cna) && !config$overwrite) {
        logwarn("ASM methylation CNA already exists. Skipping ASM methylation CNA for %s", paste0(tumor$patient_id, ":", tumor$id))
        return(asm_meth_cna)
    }

    # Get CNA solution for tumor
    cna_file <- get_fpath(tumor, config, "asm_cna")
    if (!fs::file_exists(cna_file)) {
        stop("Error. CNA file not found for tumor.")
    }
    cna <- fread_chrom(cna_file)


    # Get BAF for tumor at phased SNPs
    hap_stats <- fread_chrom(get_fpath(tumor, config, "asm_hap_stats"))

    # Get CG-hap map
    phase_map <- fread_chrom(get_fpath(tumor, config, "asm_phase_map"))

    # Get allele-specific bulk methylation
    asm_meth <- fread_chrom(get_fpath(tumor, config, "asm_meth"))

    # Overlap three datasets from ASM counter. For each CpG, return asm, BAF, cna and methylation.
    asm_hap <- merge_asm_hap(asm_meth, hap_stats, phase_map)
    # Overlap DNA methylation with CNA.
    asm_hap_cna <- overlap_meth_cna(asm_hap, cna)

    # Assign each phased CpG to a CNA state
    amc <- assign_asm_cna(asm_hap_cna)

    # Save file to system in expected location
    fs::dir_create(fs::path_dir(asm_meth_cna))
    data.table::fwrite(amc, asm_meth_cna)
    return(asm_meth_cna)
}

# TODO: Why am I not able to map all phase_map CpGs to asm_meth CpGs? It seems piledup sites are missing from the phase map (should be otherway around)

cmain_asm_deconvolve <- function(tumor, infiltrates, config) {
    # Load tumor and normal methylation
    t_meth <- fread_chrom(get_fpath(tumor, config, "asm_meth_cna"))
    n_meth <- fread_chrom(get_fpath(infiltrates, config, "asm_meth"))

    # Combine objects
    n_meth <- dplyr::rename_with(n_meth, ~ paste0(.x, "_i"), !matches("chrom|start|end"))
    setkey(n_meth, chrom, start, end)
    setkey(t_meth, chrom, start, end)
    meth_c <- merge(t_meth, n_meth, all.x = T)

    # Deconvolve ref and alt
    loginfo("Deconvolving ASM")
    meth_c <- deconvolve_asm_methylation(meth_c)

    # Filter: CN=0
    # Bulk filters not yet implemented: effective cov_t>= 3, is.na(mt-raw)
    meth_c <- meth_c[nA + nB != 0, ]

    loginfo("Calculating ASM HDI")
    # Calculate m_t HDI # parallel, long-running function
    meth_c <- calculate_asm_m_t_hdi(meth_c, config$n_cores)

    outfile <- get_fpath(tumor, config, "asm_meth_pure")
    fs::dir_create(fs::path_dir(outfile))
    data.table::fwrite(meth_c, outfile)
    return(outfile)
}

# Helper functions ----
deconvolve_asm_methylation <- function(meth_c) {
    # Deconvolve methylation for ref
    meth_c[
        ,
        ref_m_t_raw := calculate_mt(
            ref_m, ref_m_i, purity, ref_CN
        )
    ]

    # Deconvolve methylation for alt
    meth_c[
        ,
        alt_m_t_raw := calculate_mt(
            alt_m, alt_m_i, purity, alt_CN
        )
    ]

    # Correct pure tumour methylation rates set outside 0 and 1 after deconvolution
    correct_meth <- function(x) {
        data.table::fcase(
            x < 0, 0,
            x > 1, 1,
            rep(TRUE, length(x)), x
        )
    }
    meth_c[, ref_m_t := correct_meth(ref_m_t_raw)]
    meth_c[, alt_m_t := correct_meth(alt_m_t_raw)]

    # Calculate tumour coverage by deconvolution
    meth_c[, ref_cov_t := calculate_mt_cov(ref_total_counts_m, purity, ref_CN)]
    meth_c[, alt_cov_t := calculate_mt_cov(alt_total_counts_m, purity, alt_CN)]

    return(meth_c)
}

calculate_asm_m_t_hdi <- function(meth_c, n_cores, itersplit = 1e5) {
    # Split into tables of length given by itersplit for parallel processing
    inp_len <- nrow(meth_c)
    split_factor <- make_split_factor(inp_len, itersplit)
    msplit <- iterators::isplit(meth_c, split_factor)

    # Calculate HDI for both alleles
    doParallel::registerDoParallel(cores = n_cores)
    hdi <- foreach(x = msplit, .combine = "rbind") %dopar% {
        v <- x$value

        # Set empty data table to store results
        res <- data.frame(
            matrix(nrow = nrow(v), ncol = 0)
        )

        # Get table of only meth_c values eligible for HDI calculation (i.e. counts present)
        ix_asm_hdi <- sel_asm_hdi_pass(v, "ref") # Select sites eligible for HDI
        ref_hdi <- calculate_asm_hdi(v[ix_asm_hdi, ], "ref") # Calculate HDI
        res[ix_asm_hdi, colnames(ref_hdi)] <- data.frame(ref_hdi) # Assign HDI at eligible sites

        # Repeat as above for alt allele
        ix_asm_hdi <- sel_asm_hdi_pass(v, "alt") # Select sites eligible for HDI
        alt_hdi <- calculate_asm_hdi(v[ix_asm_hdi, ], "alt") # Calculate HDI
        res[ix_asm_hdi, colnames(alt_hdi)] <- data.frame(alt_hdi) # Assign HDI at eligible sites

        # Return result for binding
        return(res)
    }
    doParallel::stopImplicitCluster()

    meth_c <- cbind(meth_c, hdi)

    return(meth_c)
}

sel_asm_hdi_pass <- function(x, allele) {
    # Helper function to select sites eligible for HDI calculation for a signle allele

    # Return TRUE if allele field has counts, pure, count_i, meth_i and CN data present, otherwise return false
    # Use complete.cases function to streamline
    col_prefix <- ifelse(allele == "ref", "ref_", "alt_")
    fields <- paste0(
        col_prefix,
        c("total_counts_m", "m_t", "total_counts_m_i", "m_i", "CN")
    )
    complete.cases(x[, ..fields])
}

calculate_asm_hdi <- function(meth_c, allele) {
    # Helper function to calculate HDI for a single allele
    # Set fields based on allele
    col_prefix <- ifelse(allele == "ref", "ref_", "alt_")
    counts <- paste0(col_prefix, "total_counts_m")
    pure <- paste0(col_prefix, "m_t")
    counts_i <- paste0(col_prefix, "total_counts_m_i")
    meth_i <- paste0(col_prefix, "m_i")
    CN <- paste0(col_prefix, "CN")

    # Generate inputs for HDI calculation
    M <- round(meth_c[[counts]] * meth_c[[pure]], 0)
    UM <- meth_c[[counts]] - M
    M_n <- round(meth_c[[counts_i]] * meth_c[[meth_i]], 0)
    UM_n <- meth_c[[counts_i]] - M_n
    hdi <- vec_HDIofMCMC_mt(M, UM, M_n, UM_n, meth_c[["purity"]], meth_c[[CN]], credMass = 0.95)
    colnames(hdi) <- paste0(col_prefix, c("m_t_low", "m_t_high"))
    return(hdi)
}

cmain_asm_ss_dmps <- function(sample, config) {
    # Get pure methylation if it exists, else get bulk methylation
    asm_file <- get_fpath(sample, config, "asm_meth_pure")
    if (!file.exists(asm_file)) {
        asm_file <- get_fpath(sample, config, "asm_meth")
    }

    # Skip if output file exists
    out_file <- get_fpath(sample, config, "asm_ss_dmp")
    if (file.exists(out_file)) {
        loginfo("Skipping %s", sample$id)
        return()
    }

    # TODO: refactor function
    loginfo("Running ASM DMP calls for %s", sample$id)
    #  Calculate AS-DMP within-sample, including CAMDAC where available
    dt <- fread_chrom(asm_file)

    # Params
    # TODO: Move to config
    itersplit <- 1e5
    effect_size <- 0.2
    prob <- 0.99

    # Calculate differential methylation given effect size
    asm_b_diff <- dt[["ref_m"]] - dt[["alt_m"]]

    # Calculate DMP probability from bulk data
    # Probabilities are calculated on bulk counts
    M_alt <- round(dt[["alt_total_counts_m"]] * dt[["alt_m"]], 0)
    UM_alt <- dt[["alt_total_counts_m"]] - M_alt
    M_ref <- round(dt[["ref_total_counts_m"]] * dt[["ref_m"]], 0)
    UM_ref <- dt[["ref_total_counts_m"]] - M_ref

    if ("ref_m_t" %in% colnames(dt)) {
        asm_t_diff <- dt[["ref_m_t"]] - dt[["alt_m_t"]]
    } else {
        asm_t_diff <- NULL
    }

    # Make DMP call
    dmp_call <- dmp_call_pipe(asm_b_diff, M_ref, UM_ref, M_alt, UM_alt, asm_t_diff, effect_size, prob, itersplit)

    # Reformat and merge
    asm_names <- c("prob_DMP", "asm_b_diff", "asm_DMP_b", "asm_t_diff", "asm_DMP_t")
    names(dmp_call) <- asm_names[1:ncol(dmp_call)]
    dt <- cbind(dt, dmp_call)

    # Return data
    ss_dmp_out <- get_fpath(sample, config, "asm_ss_dmp")
    data.table::fwrite(dt, ss_dmp_out, quote = FALSE, na = "NA")
    return(ss_dmp_out)
}

cmain_asm_dmps <- function(sample, origin, config) {
    # Calculate AS-DMP between-samples
    asm_file <- get_fpath(sample, config, "asm_meth_pure")
    origin_file <- get_fpath(origin, config, "asm_meth")

    # Skip if output file exists
    out_file <- get_fpath(sample, config, "asm_dmp")
    if (file.exists(out_file)) {
        loginfo("Skipping ASM DMP calls for %s against %s", sample$id, origin$id)
        return(out_file)
    }

    loginfo("Running ASM DMP calls for %s against %s", sample$id, origin$id)

    # TODO: move to config
    effect_size <- 0.2
    prob <- 0.99
    itersplit <- 1e5

    #  Calculate AS-DMP within-sample, including CAMDAC where available
    abb <- fread_chrom(asm_file)
    ori <- fread_chrom(origin_file)
    # Rename non chrom start end names with dplyr
    ori <- dplyr::rename_with(ori, ~ paste0(.x, "_o"), .cols = !matches("chrom|start|end"))

    # Overlap datasets
    setkey(abb, chrom, start, end)
    dt <- merge(abb, ori, by = c("chrom", "start", "end"), all.x = TRUE)

    # Run DMP callin for ref allele
    mbdiff <- dt[["ref_m"]] - dt[["ref_m_o"]]
    M <- round(dt[["ref_total_counts_m_o"]] * dt[["ref_m_o"]], 0)
    UM <- dt[["ref_total_counts_m_o"]] - M
    M_n <- round(dt[["ref_total_counts_m"]] * dt[["ref_m"]], 0)
    UM_n <- dt[["ref_total_counts_m"]] - M_n
    mtdiff <- dt[["ref_m_t"]] - dt[["ref_m_o"]]

    # Make ref DMP call and merge
    loginfo("ASM: Calling tumor-normal REF DMPs")
    dmp_call <- dmp_call_pipe(mbdiff, M, UM, M_n, UM_n, mtdiff, effect_size, prob, itersplit)
    asm_names <- c("prob_ref_DMP", "ref_m_diff", "ref_DMP_b", "ref_m_t_diff", "ref_DMP_t")
    names(dmp_call) <- asm_names[1:ncol(dmp_call)]
    dt <- cbind(dt, dmp_call)

    # Run DMP calling for alt allele
    mbdiff <- dt[["alt_m"]] - dt[["alt_m_o"]]
    M <- round(dt[["alt_total_counts_m_o"]] * dt[["alt_m_o"]], 0)
    UM <- ori[["alt_total_counts_m"]] - M
    M_n <- round(dt[["alt_total_counts_m"]] * dt[["alt_m"]], 0)
    UM_n <- dt[["alt_total_counts_m"]] - M_n
    mtdiff <- dt[["alt_m_t"]] - dt[["alt_m_o"]]

    # Make alt DMP calla nd merge
    loginfo("ASM: Calling tumor-normal ALT DMPs")
    dmp_call <- dmp_call_pipe(mbdiff, M, UM, M_n, UM_n, mtdiff, effect_size, prob, itersplit)
    asm_names <- c("prob_alt_DMP", "alt_m_diff", "alt_DMP_b", "alt_m_t_diff", "alt_DMP_t")
    names(dmp_call) <- asm_names[1:ncol(dmp_call)]
    dt <- cbind(dt, dmp_call)

    # Save data
    dmp_out <- get_fpath(sample, config, "asm_dmp")
    data.table::fwrite(dt, dmp_out, quote = FALSE, na = "NA")
    return(dmp_out)
}
