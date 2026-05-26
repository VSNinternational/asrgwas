library(asreml)
asreml.options(trace = FALSE)
suppressWarnings(gwas.data <<- pre.gwas(
  pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
  geno.data = geno.apricot, map.data = map.apricot, maf = 0.05,
  marker.callrate = 0.40, impute = TRUE))
Kinv_wrong <- gwas.data$Kinv
attributes(Kinv_wrong[[1]])$colNames[1] <- "wrong name"
Q_wrong <- gwas.data$Q
rownames(Q_wrong)[1] <- "wrong name"
####################################
test_that("match.data simple call", {
  expect_no_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      message = TRUE
      )
    )
  expect_no_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data, gwas.data$geno.data),
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      message = TRUE
    )
  )
})
####################################
test_that("test traps on match.data", {
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data[-1,],
      geno.data = list(gwas.data$geno.data),
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      message = TRUE
    )
  )
  expect_error(
    match.data(
      pheno.data = "",
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = "",
      geno.data = list(gwas.data$geno.data),
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(""),
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = "",
      Q = gwas.data$Q,
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = list(gwas.data$Kinv, gwas.data$Kinv),
      Q = gwas.data$Q,
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = list(Kinv_wrong),
      Q = gwas.data$Q,
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = list(gwas.data$Kinv, gwas.data$Kinv),
      Q = "",
      message = TRUE)
  )
  expect_error(
    match.data(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      map.data = gwas.data$map.data,
      geno.data = list(gwas.data$geno.data),
      Kinv = list(gwas.data$Kinv),
      Q = Q_wrong,
      message = TRUE)
  )
})
