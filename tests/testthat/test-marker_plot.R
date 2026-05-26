library(asreml)
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
####################################
test_that("marker.plot simple call", {
  expect_no_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = c("Sucrose", "Ethylene"),
      geno.data.sel = geno.data.sel,
      type.plot = "boxplot",
      family = "gaussian",
      sample.label = TRUE,
      sample.width = TRUE,
      maf.colour = TRUE,
      facet.dim = c(2,4),
      scales = NULL,
      message = TRUE)
  )
  expect_no_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      geno.data.sel = geno.data.sel,
      type.plot = "violin",
      family = "gaussian",
      sample.label = FALSE,
      sample.width = TRUE,
      maf.colour = TRUE,
      facet.dim = NULL,
      scales = "free_y",
      message = TRUE)
  )
  expect_no_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = c("Sucrose", "Ethylene"),
      geno.data.sel = geno.data.sel,
      type.plot = "violin",
      family = "binomial",
      sample.label = TRUE,
      sample.width = TRUE,
      maf.colour = TRUE,
      facet.dim = c(4,2),
      scales = "free",
      message = TRUE)
  )
  expect_no_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Ethylene",
      geno.data.sel = geno.data.sel,
      type.plot = "boxplot",
      family = "binomial",
      sample.label = TRUE,
      sample.width = TRUE,
      maf.colour = TRUE,
      facet.dim = NULL,
      scales = "free_x",
      message = TRUE)
  )
})
####################################
test_that("test traps on marker.plot", {
  expect_error(
    marker.plot(
      pheno.data = "", indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      sample.label = TRUE, sample.width = TRUE, maf.colour = TRUE,
      scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      sample.label = TRUE, sample.width = TRUE, maf.colour = TRUE,
      scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = "", geno.data.sel = geno.data.sel,
      sample.label = TRUE, sample.width = TRUE, maf.colour = TRUE,
      scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = "",
      sample.label = TRUE, sample.width = TRUE, maf.colour = TRUE,
      scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      type.plot = "", sample.label = TRUE, sample.width = TRUE,
      maf.colour = TRUE, scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      family = "", sample.label = TRUE, sample.width = TRUE,
      maf.colour = TRUE, scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      sample.label = "", sample.width = TRUE, maf.colour = TRUE,
      scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      sample.label = TRUE, sample.width = "", maf.colour = TRUE,
      scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      sample.label = TRUE, sample.width = TRUE, maf.colour = "",
      scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      facet.dim = "", sample.label = TRUE, sample.width = TRUE,
      maf.colour = TRUE, scales = "free")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      sample.label = TRUE, sample.width = TRUE, maf.colour = TRUE,
      scales = "")
    )
  expect_error(
    marker.plot(
      pheno.data = gwas.data$pheno.data, indiv = "Ind",
      resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
      sample.label = TRUE, sample.width = TRUE, maf.colour = TRUE,
      scales = "free", message = "")
    )
})
