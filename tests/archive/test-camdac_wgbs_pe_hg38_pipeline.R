# context("Test wrapper functions for CAMDAC WGBS")
# 
# test_that("Alleles are counted", {
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
#   allele_counts_file <- cmain_count_alleles(tumour_obj, run_config)
#   allele_counts_file <- cmain_count_alleles(normal_obj, run_config)
#   
#   expect_equal(
#     allele_counts_file,
#     "./runs/pgp-test/Allelecounts/pgp-001/pgp-test.pgp-001.SNPs.CpGs.all.sorted.csv.gz"
#   )
#   
#   acf = readRDS(allele_counts_file)
#   expect_equal(
#     names(acf),
#     c("CHR","chrom","start","end","width","POS","ref","alt","alt_counts","ref_counts",
#     "total_counts","BAF","total_depth","other_counts","all_counts","M","UM",
#     "total_counts_m","m","Af","Ar","Cf","Cr","Tf","Tr","Gf","Gr","CAr","TGf",
#     "CGf","CGr","CCGG","mq")
#   )
#   rm(acf)
# })
# 
