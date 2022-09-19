# # Test methylation processing functions for pe hg38 wgbs data
# ##load_all(); test_file("tests/testthat/test-camdac_wgbs_methylation.R")
# 
# # Setup methylation test data
# tumour_obj <- create_camdac_sample(
#   "pgp-test", "XY", "pgp-001", "tumour", "./data/pgp_3255_101_cpgs.bam"
# )
# 
# normal_obj <- create_camdac_sample(
#   "pgp_3255", "XY", "pgp-002", "normal", "./data/pgp_3255_103_cpgs.bam"
# )
# 
# run_config <- create_camdac_config(
#   camdac_refs="./data/pipeline_files", outdir='./runs/',
#   build="hg38", bsseq="wgbs", bsseq_lib="pe", n_cores=2
# )
# 
# ac_file <- "data/pgp-test.pgp-001.SNPs.CpGs.all.sorted.rds"
# 
# test_that("some_function", {
# 
#   tumour_obj <- create_camdac_sample(
#     patient_id="pgp-test",
#     patient_sex="XY",
#     sample_id="pgp-000",
#     sample_type="tumour",
#     bam_file="./data/bam/pgp_3255_101_cpgs.bam"
#   )
#   
#   normal_obj <- create_camdac_sample(
#     patient_id="pgp-test",
#     patient_sex="XY",
#     sample_id="pgp-001",
#     sample_type="normal",
#     bam_file="./data/bam/pgp_3255_103_cpgs.bam"
#   )  
#   
#   run_config <- create_camdac_config(
#     camdac_refs="./data/pipeline_files",
#     outdir="./runs/",
#     build="hg38",
#     bsseq="wgbs",
#     bsseq_lib="pe",
#     n_cores=5
#   )
#   
#   tumour_ac <- build_ouput
#   expect_true(fs::file_exists(build_output))
#     
# }
# )
# 
# 
# # test_that("tumour and normal ac files exist", {
# #   tumour_ac_file <- build_output_name(tumour_sample, test_meth_config, "allele_counts")
# #   if (!fs::file_exists(tumour_ac_file)){
# #     allele_counts_file <- cmain_count_alleles(tumour_sample, test_meth_config)
# #   }
# #   
# #   normal_ac_file <- build_output_name(normal_sample, test_meth_config, "allele_counts")
# #   if (!fs::file_exists(normal_ac_file)){
# #     allele_counts_file <- cmain_count_alleles(normal_sample, test_meth_config)
# #   }
# #    
# #   expect_true(
# #     all(fs::file_exists(tumour_ac_file), fs::file_exists(tumour_ac_file))
# #   )
# # })
# 
# # test_that("allele counts formatted for methylation analysis", {
# #   normal_ac_file <- build_output_name(normal_sample, test_meth_config, "allele_counts")
# #   normal_ac <- readRDS(normal_ac_file)
# #   normal_m <- format_counts_to_methyl(normal_ac, min_cov=3)
# #   
# #   # Conditions:
# #   # No SNP sites present, methylation rates are not NA, total_counts_m above threshold
# #   #   expected new columns are present, no duplicate loci, 
# #   # 
# #   expect_true(all(!is.na(tumour_m$POS)))
# #   expect_true(all(!is.na(tumour_m$total_counts_m)))
# #   expect_equal(names(tumour_m), c("seqnames","start","end","width","strand","M_b","UM_b","m_b","cov_b","SNP_b"))
# # })
# # 
# # 
# # test_that("allele counts formatted for methylation analysis", {
# #   tumour_ac_file <- build_output_name(tumour_sample, test_meth_config, "allele_counts")
# #   tumour_ac <- readRDS(tumour_ac_file)
# #   tumour_m <- format_counts_to_methyl(tumour_ac, min_cov=3)
# #   
# #   # Conditions:
# #   # No SNP sites present, methylation rates are not NA, total_counts_m above threshold
# #   #   expected new columns are present, no duplicate loci, 
# #   # 
# #   expect_true(all(!is.na(tumour_m$POS)))
# #   expect_true(all(!is.na(tumour_m$total_counts_m)))
# #   expect_equal(names(tumour_m), c("seqnames","start","end","width","strand","M_b","UM_b","m_b","cov_b","SNP_b"))
# # })
# # 
