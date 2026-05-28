testthat::skip_if_not_installed("asreml")
suppressPackageStartupMessages(library(asreml))
asreml.options(trace = FALSE)
gwas_data <- ASRgwas::map.apricot
gwas_data_chr1 <- gwas_data[gwas_data$chrom == "Pp01",]
############################################
test_that("marker.position simple call", {
  expect_no_error(
    marker.position(map.data = gwas_data_chr1,
                    chrom = "chrom",
                    pos = "pos",
                    collate = TRUE,
                    padding = 0,
                    message = TRUE)
    )
})
############################################
