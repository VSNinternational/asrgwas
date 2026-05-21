test_that("call.generator simple call", {
  expect_no_error(
    ASRgwas:::call.generator(
      resp = "resp",
      fixedf = c("year", "block"),
      group.cols = c(10,12),
      cov = "stand",
      randomf = c("gen:year"),
      randomf.vm = "gen",
      vm.matrix.name = "Kinverse",
      data.name = "phenotypic.data",
      residual = "year",
      weights = "weights",
      family = "gaussian",
      workspace = "4Gb"
      )
  )
  expect_no_error(
    ASRgwas:::call.generator(
      resp = "resp",
      fixedf = c("year", "block"),
      group.cols = c(10,12),
      cov = "stand",
      randomf = c("gen:year"),
      randomf.vm = "gen",
      vm.matrix.name = "Kinverse",
      data.name = "phenotypic.data",
      residual = NULL,
      weights = NULL,
      family = "gaussian",
      workspace = "4Gb"
    )
  )
  expect_no_error(
    ASRgwas:::call.generator(
      resp = "resp",
      fixedf = c("year", "block"),
      group.cols = c(10,12),
      cov = "stand",
      randomf = c("gen:year"),
      randomf.vm = "gen",
      vm.matrix.name = "Kinverse",
      data.name = "phenotypic.data",
      residual = NULL,
      weights = NULL,
      family = "binomial",
      dispersion = 1,
      total = NULL,
      workspace = "4Gb"
    )
  )
  expect_no_error(
    ASRgwas:::call.generator(
      resp = "resp",
      fixedf = c("year", "block"),
      group.cols = c(10,12),
      cov = "stand",
      randomf = c("gen:year"),
      randomf.vm = NULL,
      vm.matrix.name = "Kinverse",
      data.name = "phenotypic.data",
      residual = NULL,
      weights = NULL,
      family = "binomial",
      dispersion = 1,
      total = "total",
      workspace = "4Gb"
    )
  )
  expect_no_error(
    ASRgwas:::call.generator(
      resp = "resp",
      fixedf = c("year", "block"),
      group.cols = c(10,12),
      cov = "stand",
      randomf = NULL,
      randomf.vm = NULL,
      vm.matrix.name = "Kinverse",
      data.name = "phenotypic.data",
      residual = NULL,
      weights = NULL,
      family = "binomial",
      dispersion = 1,
      total = "total",
      workspace = "4Gb"
    )
  )
})
