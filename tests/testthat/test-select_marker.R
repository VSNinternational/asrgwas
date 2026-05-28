testthat::skip_if_not_installed("asreml")
suppressPackageStartupMessages(library(asreml))
asreml.options(trace = FALSE)
suppressWarnings(gwas.data <<- pre.gwas(
  pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
  geno.data = geno.apricot, map.data = map.apricot, maf = 0.05,
  marker.callrate = 0.40, impute = TRUE))
gwasA <<- gwas.asreml(
  pheno.data = gwas.data$pheno.data, resp = "Sucrose", gen = "Ind",
  fixedf = "Lots", residual = "Lots",
  Kinv = gwas.data$Kinv, Q = gwas.data$Q, npc = 3,
  geno.data = gwas.data$geno.data, map.data = gwas.data$map.data,
  pvalue.thr = 0.0005, bonferroni = FALSE)
sign.markers <- gwasA$gwas.sel$marker
if (is.null(sign.markers) || length(sign.markers) == 0) {
  sign.markers <- head(gwasA$gwas.all$marker, 5)
}
geno.data.sel <- as.matrix(gwas.data$geno.data[, sign.markers, drop = FALSE])
gwasA$mod$vparameters
sign.markers2 <- gwasA$gwas.all[1:5,]$marker
geno.data.sel2 <- as.matrix(gwas.data$geno.data[, sign.markers2, drop = FALSE])
####################################
test_that("select.marker simple call", {
  expect_no_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = geno.data.sel,
      ref.vc = 1,
      pvalue.thr = 0.1,
      maxiter.update = 5,
      workspace = "1Gb",
      message = TRUE
      )
  )
  expect_no_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = geno.data.sel2,
      ref.vc = 1,
      pvalue.thr = 0.1,
      maxiter.update = 5,
      workspace = "1Gb",
      message = TRUE
    )
  )
})
####################################
test_that("test traps on select.marker", {
  expect_error(
    select.marker(
      gwas.object = "",
      geno.data.sel = geno.data.sel,
      ref.vc = 1,
      pvalue.thr = 0.1,
      maxiter.update = 5,
      workspace = "1Gb",
      message = TRUE
    )
  )
  expect_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = "",
      ref.vc = 1,
      pvalue.thr = 0.1,
      maxiter.update = 5,
      workspace = "1Gb",
      message = TRUE
    )
  )
  expect_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = geno.data.sel,
      ref.vc = "",
      pvalue.thr = 0.1,
      maxiter.update = 5,
      workspace = "1Gb",
      message = TRUE
    )
  )
  expect_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = geno.data.sel,
      ref.vc = 1,
      pvalue.thr = "",
      maxiter.update = 5,
      workspace = "1Gb",
      message = TRUE
    )
  )
  expect_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = geno.data.sel,
      ref.vc = 1,
      pvalue.thr = 0.1,
      maxiter.update = "",
      workspace = "1Gb",
      message = TRUE
    )
  )
  expect_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = geno.data.sel,
      ref.vc = 1,
      pvalue.thr = 0.1,
      maxiter.update = 5,
      workspace = "",
      message = TRUE
    )
  )
  expect_error(
    select.marker(
      gwas.object = gwasA,
      geno.data.sel = geno.data.sel,
      ref.vc = 1,
      pvalue.thr = 0.1,
      maxiter.update = 5,
      workspace = "1Gb",
      message = ""
    )
  )
})
