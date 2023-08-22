test_that("CAMDAC runs with battenberg in hg19 mode", {
    # CNA caller test config
    config_hg19 <- CamConfig(
        outdir = "./result_hg19_pe",
        bsseq = "wgbs",
        build = "hg19",
        lib = "pe",
        # hg19 Battenberg will fail as not enough SNPs on each hg19/hg38 chrom for haplotyping to work.
        # Hence, can't use tumor-normal and must SNP-inject?
        cna_caller = "battenberg",
        n_cores = 10,
        min_cov = 1 # Required to capture sufficient SNPs from test
    )
    config_hg19$regions = system.file("testdata", "test_wgbs_segments.bed", package = "CAMDAC")

    # Mock het SNP selection. As we're not using true tumor-normal pairs,
    # normal SNP selection on BAF 0.2 and 0.8 will fail to yield sufficient SNPs
    local_mocked_bindings(
        select_heterozygous_snps = function(tsnps, ...){
            return(tsnps)
        }
    )

    local_mocked_bindings(
        get.chrom.names = function(...){
            c("3", "9")
        },
        .package="Battenberg"
    )

    # Try again, simply replacing tsnps file
    tsnps_file <- system.file("testdata", "test_tsnps.csv.gz", package = "CAMDAC")
    attach_output(tumor, config_hg19, "tsnps", tsnps_file)

    stdout <- testthat::capture_output(
        pipeline(tumor, germline = normal, infiltrates = normal, origin = normal, config_hg19)
    )

    expr <- "pipeline complete"

    testthat::expect_true(
        stringr::str_detect(stdout, expr)
    )

})
