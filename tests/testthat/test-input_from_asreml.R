library(asreml)
asreml.options(trace = FALSE)
model1 <<- asreml::asreml(yield ~ Variety*Nitrogen,
                          random = ~ Blocks/Wplots,
                          residual = ~units,
                          data = oats)
model2 <<- asreml::asreml(yield ~ Variety*Nitrogen,
                          random = ~ Blocks/Wplots,
                          residual = ~dsum(~ units | Blocks),
                          data = oats)
####################################
test_that("input.from.asreml simple call", {
  expect_no_error(
    input.from.asreml(mod = model1, parent.frame = 1)
  )
  expect_no_error(
    input.from.asreml(mod = model2, parent.frame = 1)
  )
})
