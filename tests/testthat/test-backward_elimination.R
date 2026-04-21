library(asreml)
asreml.options(trace = FALSE)
model <<- asreml::asreml(yield ~ Variety*Nitrogen,
                        random = ~ Blocks/Wplots,
                        na.action = list(y='omit', x='omit'),
                        data = oats)
####################################
test_that("backward.elimination simple call", {
  expect_no_error(
    backward.elimination(
      mod = model,
      data = oats,
      try.variables = c("Row", "Column"),
      pvalue.thr = 0.1,
      maxiter.update = 10
      )
  )
})
####################################
test_that("test traps on backward.elimination", {
  expect_error(
    backward.elimination(
      mod = "",
      data = oats,
      try.variables = c("Row", "Column"),
      pvalue.thr = 0.1,
      maxiter.update = 10
    )
  )
  expect_error(
    backward.elimination(
      mod = model,
      data = "",
      try.variables = c("Row", "Column"),
      pvalue.thr = 0.1,
      maxiter.update = 10
    )
  )
  expect_error(
    backward.elimination(
      mod = model,
      data = oats,
      try.variables = "",
      pvalue.thr = 0.1,
      maxiter.update = 10
    )
  )
  expect_error(
    backward.elimination(
      mod = model,
      data = oats,
      try.variables = c("Row", "Column"),
      pvalue.thr = "",
      maxiter.update = 10
    )
  )
  expect_error(
    backward.elimination(
      mod = model,
      data = oats,
      try.variables = c("Row", "Column"),
      pvalue.thr = 0.1,
      maxiter.update = ""
    )
  )
})
