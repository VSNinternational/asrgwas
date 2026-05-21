library(asreml)
asreml.options(trace = FALSE)
model <<- asreml(fixed=Sucrose~1+Ind,
                 random = ~Lots,
                 residual=~idv(units),
                 data=ASRgwas::pheno.apricot)
############################################
test_that("v_aux simple call", {
  expect_no_error(
    check.special("variance", "Lots", model$G.param, "id")
  )
  expect_no_error(
    check.special("Lots", "Lots", model$G.param, "id")
  )
  expect_no_error(
    g.asreml("variance",
             "Lots",
             model$G.param,
             ASRgwas:::check.special("variance", "Lots", model$G.param, "id"),
             cond.fac = "")
    )
})
############################################
test_that("test traps on v_aux", {
  expect_error(
    check.special("Lots", "Lots", model$G.param, "",residual = TRUE)
  )
})
