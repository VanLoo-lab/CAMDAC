test_that("CNA calls can be made without a germline sample", {
    # CNA caller test config
    config_c <- CamConfig(
        outdir = "./result_to",
        bsseq = "wgbs",
        build = "hg38",
        lib = "pe",
        regions = regions,
        n_cores = 30,
        cna_caller = "ascat",
        min_cov = 1 # Required to capture sufficient SNPs from test
    )
    #withr::defer(fs::dir_delete(config_c$outdir))

    # Setup germline versions: no data, SNP pos only, SNP pos with counts
    germline <- normal
    germline_to <- NULL 

    germline_pos <- CamSample(id = "G", sex = "XY")
    germline_pos_file <- system.file("testdata", "test.to.norm_pos.csv.gz", package = "CAMDAC")
    attach_output(germline_pos, config_c, "snps", germline_pos_file)

    germline_count <- CamSample(id= "GP", sex= "XY")
    germline_count_file <- system.file("testdata", "test.to.norm_pos_count.csv.gz", package = "CAMDAC")
    attach_output(germline_count, config_c, "snps", germline_count_file)

    # Run allele counting
    preprocess(list(
        tumor, germline, germline_to, germline_count, germline_pos
    ), config_c)

    # Mock gc annotation
    local_mocked_bindings(
        select_heterozygous_snps = function(tsnps, ...){
            tsnps
        }
    )

    # Set output names expectation
    exp_names = c("chrom", "POS", "total_counts", "total_depth", "ref", "alt", "BAF", "BAFr",
     "total_depth_n", "total_counts_n", "BAF_n", "BAFr_n", "LogR_n")

    # Tumor only
    tsnps_to = cmain_bind_snps(tumor, germline_to, config_c)
    fs::file_exists(tsnps_to)
    fs::file_delete(c(tsnps_to))

    # SNP only
    tsnps_pos = cmain_bind_snps(tumor, germline_pos, config_c)
    fs::file_exists(tsnps_pos)
    fs::file_delete(c(tsnps_pos))
    
    # SNP with counts
    tsnps_count = cmain_bind_snps(tumor, germline_count, config_c)
    fs::file_exists(tsnps_count)
    fs::file_delete(c(tsnps_count))

})
