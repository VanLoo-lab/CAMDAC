test_that("ASM pipeline runs", {
  b_tumor <- system.file("testdata", "tumor.bam", package = "CAMDAC")
  b_normal <- system.file("testdata", "normal.bam", package = "CAMDAC")

  tumor <- CamSample(id="T", sex="XY", bam=b_tumor)
  normal <- CamSample(id="N", sex="XY", bam=b_normal)
  config <- CamConfig(
    outdir="./result_asm_full", bsseq="wgbs", lib="pe",
    build="hg38", n_cores=3, min_cov=1, cna_caller="ascat")

  asm_pipeline(
      tumor=tumor,
      germline=normal,
      infiltrates=normal,
      origin=normal,
      config=config
  )
})
