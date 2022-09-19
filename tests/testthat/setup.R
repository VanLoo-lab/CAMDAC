
# Global objects that may be used by tests.
# TODO: Determine whether these fixtures are required
# tumor <- create_camdac_sample(
#   patient_id="TEST",
#   patient_sex="XX",
#   sample_id="T",
#   sample_type="tumour",
#   bam_file=system.file("extdata", "NA18939_test_v3.bam", package = "CAMDAC")
# )
# 
# normal <- create_camdac_sample(
#   patient_id="TEST",
#   patient_sex="XX",
#   sample_id="N",
#   sample_type="tumour",
#   bam_file="/camp/home/mensahn/pvs/camdac_pe_test/g1000_test/NA20502_test_v3.bam"
# )
# 
# config <- create_camdac_config(
#   camdac_refs="/camp/project/proj-vanloo/secure/working/nmensah/pipeline_files",
#   outdir="/camp/project/proj-vanloo/secure/working/nmensah/camdac_pe_test/g1000_test/test_result_full",
#   build="hg38",
#   bsseq="wgbs",
#   bsseq_lib="pe",
#   n_cores=8,
#   n_seg_split=50,
#   min_cov = 1,
#   ascat_rho_manual = 0.8, # Required for test data battenberg else no solution found
#   ascat_psi_manual = 2.8 # Required for test data battenberg else no solution found
# )

# Cleanup, as presented in https://testthat.r-lib.org/articles/test-fixtures.html
withr::defer(teardown_env())