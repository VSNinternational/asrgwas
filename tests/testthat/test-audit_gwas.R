suppressWarnings(gwas.data <<- pre.gwas(
  pheno.data = pheno.apricot, indiv = "Ind",
  resp = "Sucrose", geno.data = geno.apricot))
####################################
test_that("audit.gwas simple call", {
  expect_no_error(
    audit.gwas(
      pheno.data = gwas.data$pheno.data,
      indiv = "Ind",
      resp = "Sucrose",
      Kinv = gwas.data$Kinv,
      geno.data = gwas.data$geno.data,
      map.data = NULL,
      Q = gwas.data$Q,
      message = TRUE
      )
  )
})
#####################################
test_that("test traps on audit.gwas", {
  expect_error(
    audit.gwas(
      pheno.data = "", resp = "Sucrose", indiv = "Ind",
      Kinv = data$Kinv, geno.data = data$geno.data, Q = data$Q)
    )
  expect_error(
    audit.gwas(
      pheno.data = data$pheno.data, resp = "", indiv = "Ind",
      Kinv = data$Kinv, geno.data = data$geno.data, Q = data$Q)
    )
  expect_error(
    audit.gwas(
      pheno.data = data$pheno.data, resp = "Sucrose", indiv = "",
      Kinv = data$Kinv, geno.data = data$geno.data, Q = data$Q)
    )
  expect_error(
    audit.gwas(
      pheno.data = data$pheno.data, resp = "Sucrose", indiv = "Ind",
      Kinv = "", geno.data = data$geno.data, Q = data$Q)
    )
  expect_error(
    audit.gwas(
      pheno.data = data$pheno.data, resp = "Sucrose", indiv = "Ind",
      Kinv = data$Kinv, geno.data = "", Q = data$Q)
    )
  expect_error(
    audit.gwas(
      pheno.data = data$pheno.data, resp = "Sucrose", indiv = "Ind",
      Kinv = data$Kinv, geno.data = data$geno.data, Q = "")
    )
})
###################################
  ## Expect "audit.gwas" type.
