test_that("pre.gwas simple call", {
  expect_warning(
    pre.gwas(
      pheno.data = pheno.apricot,
      indiv = "Ind",
      resp = "Sucrose",
      geno.data = geno.apricot)
  )
  expect_no_error(
    pre.gwas(
      pheno.data = NULL,
      indiv = NULL,
      geno.data = NULL,
      resp = NULL,
      map.data = NULL,
      )
  )
})
