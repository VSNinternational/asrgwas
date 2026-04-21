library(asreml)
asreml.options(trace = FALSE)
model <<- asreml(fixed=Sucrose~1+Ind,
                 residual=~idv(units),
                 na.action=list(x='include',y='include'),
                 data=ASRgwas::pheno.apricot)
############################################
test_that("test traps on repeatability.pev", {
  expect_error(
    repeatability.pev(
      mod = lm(yield ~ Variety*Nitrogen, data = oats),
      component = NULL
    )
  )
  expect_error(
    repeatability.pev(
      mod = NULL,
      component = NULL
    )
  )
  expect_error(
    repeatability.pev(
      mod = model,
      component = NULL
    )
  )
})
