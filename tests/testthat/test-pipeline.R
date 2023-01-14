test_that("Pipeline tumour-normal completes", {
    stdout <- testthat::capture_output(
        pipeline(tumor, germline = normal, infiltrates = normal, origin = normal, config)
    )

    expr <- "pipeline complete"

    testthat::expect_true(
        stringr::str_detect(stdout, expr)
    )
})
