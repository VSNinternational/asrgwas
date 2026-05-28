testthat::skip_if_not_installed("asreml")
suppressPackageStartupMessages(library(asreml))
asreml.options(trace = FALSE)
suppressWarnings(gwas.data <<- pre.gwas(
  pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
  geno.data = geno.apricot, map.data = map.apricot,
  maf = 0.05, heterozygosity = 0.9, Fis = 0.5,
  marker.callrate = 0.40, impute = TRUE))
caller <<- "gwas"
null_rnames <- gwas.data$geno.data
null_cnames <- gwas.data$geno.data
rownames(null_rnames) <- NULL
colnames(null_cnames) <- NULL
pheno.data2 <- gwas.data$pheno.data[,-1]
pheno.data2$Ind <- gwas.data$pheno.data$Ind
pheno.data2$Sucrose[pheno.data2$Sucrose %in% sample(pheno.data2$Sucrose, 200)] <- NA
caller <<- "gwas"
######################################
test_that("pre.generic simple call", {
  expect_no_error(
    pre.generic(
      pheno.data = NULL,
      indiv = NULL,
      resp = "Sucrose",
      weights = NULL,
      geno.data = NULL,
      map.data = NULL
    )
  )
  expect_warning(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = as.matrix(gwas.data$K$Ind),
      Q.method = "none",
      return = "pheno.data",
      message = TRUE,
      ellipsis = list(clean.diagonal = TRUE)
    )
  )
  expect_warning(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = as.matrix(gwas.data$K$Ind),
      Q.method = "K",
      return = "Kinv",
      message = TRUE,
      ellipsis = list(clean.diagonal = TRUE)
    )
  )
  expect_warning(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = as.matrix(gwas.data$K$Ind),
      Q.method = "geno.data",
      return = "pheno.data",
      message = TRUE,
      ellipsis = list(clean.diagonal = TRUE)
    )
  )
  expect_no_error(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = NULL,
      Q.method = "none",
      return = "pheno.data",
      message = TRUE,
      ellipsis = list(pruning.thr = TRUE)
    )
  )
})
########################################
test_that("test traps on pre.generic", {
  expect_error(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = as.matrix(gwas.data$Kinv$Ind),
      Q.method = "none",
      return = "pheno.data",
      message = TRUE,
      ellipsis = NULL
    )
  )
  expect_error(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = NULL,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = NULL,
      Q.method = "none",
      return = "geno.data",
      message = TRUE,
      ellipsis = NULL
    )
  )
  expect_error(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = null_rnames,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = as.matrix(gwas.data$K$Ind),
      Q.method = "none",
      return = "pheno.data",
      message = TRUE,
      ellipsis = list(clean.diagonal = TRUE)
    )
  )
  expect_error(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = null_cnames,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = as.matrix(gwas.data$K$Ind),
      Q.method = "none",
      return = "pheno.data",
      message = TRUE,
      ellipsis = list(clean.diagonal = TRUE)
    )
  )
  expect_error(
    pre.generic(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data[,-1],
      rename.markers = TRUE,
      K = as.matrix(gwas.data$K$Ind),
      Q.method = "none",
      return = "pheno.data",
      message = TRUE,
      ellipsis = NULL
      )
  )
  expect_error(
    pre.generic(
      pheno.data = pheno.data2,
      indiv = "Ind",
      resp = "Sucrose",
      weights = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      rename.markers = TRUE,
      K = as.matrix(gwas.data$K$Ind),
      Q.method = "none",
      return = "pheno.data",
      message = TRUE,
      ellipsis = NULL
    )
  )
})
