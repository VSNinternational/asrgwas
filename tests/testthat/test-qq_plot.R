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
n_rows_gwasA <- nrow(gwasA$gwas.all)
gwasA$gwas.all2 <- gwasA$gwas.all
gwasA$gwas.all2$index <- c(rep(1, floor(n_rows_gwasA/2)), rep(2, ceiling(n_rows_gwasA/2)))
gwasA$gwas.all2[2, "marker"] <- gwasA$gwas.all[1, "marker"]
####################################
test_that("qq.plot simple call", {
  expect_no_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
  )
  expect_no_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "chrom",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    qq.plot(
      gwas.table = gwasA$gwas.all2,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    qq.plot(
      gwas.table = gwasA$gwas.all2,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = "index",
      facet.main = "index",
      facet.secondary = NULL,
      facet.dim = c(1, 2),
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    qq.plot(
      gwas.table = gwasA$gwas.all2,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = "index",
      facet.main = "index",
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
    )
  )
})
####################################
test_that("test traps on qq.plot", {
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all[,3:4],
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = FALSE,
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = "",
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = "",
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = FALSE,
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = "",
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = FALSE,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = 1,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = 1,
      facet.dim = NULL,
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = "",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_warning(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = c(10, 1),
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    qq.plot(
      gwas.table = gwasA$gwas.all,
      point.colour = "green",
      point.size = 2,
      point.alpha = 0.5,
      trend.colour = "grey",
      coord.equal = TRUE,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      legend.position = "none",
      message = ""
      )
    )
})
