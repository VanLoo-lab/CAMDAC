
#' Call CAMDAC for a tumor and patient-matched normal sample
#'
#' @param patient_id character. Patient identifier
#' @param tumor_id character. Tumor sample identifier
#' @param normal_id character. Normal sample identifier
#' @param tumor_bam character. Full path to tumor bam file
#' @param normal_bam character. Full path to normal bam file
#' @param sex character. Patient sex: "XX" for female or "XY" for male
#' @param path character. Full path to CAMDAC output directory
#' @param pipeline_files character. Full path to parent directory containing CAMDAC pipeline_files
#' @param build character. Genome build: "hg19" or "hg38"
#' @param min_tumor integer. Minimum read filter for tumor samples
#' @param min_normal integer. Minimum read filter for normal samples
#' @param n_cores integer. Number of cores to use for parallel processing
#' @param mq integer. Minimum mapping quality filter
#' @export
pipeline_tumor_normal <- function(patient_id, tumor_id, normal_id, tumor_bam, normal_bam, sex, path,
                                  pipeline_files, build, min_tumor = 3, min_normal = 10,
                                  n_cores = 1, mq = 0, paired_end) {
    # Preprocess tumor and normal sample
    preprocess_sample(
        patient_id, normal_id, normal_id, normal_bam, min_tumor,
        min_normal, mq, sex, path, pipeline_files, build, n_cores, paired_end
    )
    preprocess_sample(
        patient_id, tumor_id, normal_id, tumor_bam, min_tumor,
        min_normal, mq, sex, path, pipeline_files, build, n_cores, paired_end
    )

    # Get purified methylation rate
    get_pure_tumour_methylation(
        patient_id = patient_id, sample_id = tumor_id, sex = sex,
        normal_infiltrates_proxy_id = normal_id,
        path, pipeline_files, build,
        n_cores, reseg = FALSE
    )

    # Get DMP and DMR calls
    get_differential_methylation(
        patient_id = patient_id, sample_id = tumor_id, sex = sex,
        normal_origin_proxy_id = normal_id,
        path, pipeline_files, build,
        effect_size = 0.2, prob = 0.99,
        min_DMP_counts_in_DMR = 5, min_consec_DMP_in_DMR = 4,
        n_cores, reseg = FALSE, bulk = FALSE
    )
}


preprocess_sample <- function(patient_id, sample_id, normal_id, bam_file, min_tumor,
                              min_normal, mq, sex, path, pipeline_files, build, n_cores) {
    # Run allele counter for normal sample
    for (a in 1:25) {
        get_allele_counts(
            i = a, patient_id = patient_id, sample_id = sample_id, sex = sex, bam_file, mq = mq,
            path, pipeline_files, build, n_cores, test = FALSE, paired_end
        )
    }

    # Merge allele counts
    # Set normal status based on whether sample and normal ID match
    is_normal <- ifelse(sample_id == normal_id, TRUE, FALSE)
    format_output(
        patient_id, sample_id, sex, is_normal, path, pipeline_files, build
    )

    # Create SNP files (normal) or run ASCAT (tumor)
    run_ASCAT.m(
        patient_id, sample_id, sex,
        patient_matched_normal_id = normal_id,
        path, pipeline_files, build,
        min_normal, min_tumor,
        n_cores, reference_panel_coverage = NULL
    )

    # Process methylation info for copy number profiling and plot summary.
    run_methylation_data_processing(
        patient_id, sample_id,
        normal_infiltrates_proxy_id = normal_id,
        normal_origin_proxy_id = normal_id,
        path, min_normal, min_tumor, n_cores,
        reference_panel_normal_infiltrates = NULL,
        reference_panel_normal_origin = NULL
    )
}
