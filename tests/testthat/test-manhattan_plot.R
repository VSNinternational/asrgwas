testthat::skip_if_not_installed("asreml")
suppressPackageStartupMessages(library(asreml))
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
tag.table <- head(gwasA$gwas.all[order(gwasA$gwas.all$p.value), ], 5)
tag.table$tag <- paste("Putative gene", seq_len(nrow(tag.table)))
n_rows_gwasA <- nrow(gwasA$gwas.all)
gwasA$gwas.all2 <- gwasA$gwas.all
gwasA$gwas.all2$index <- c(rep(1, floor(n_rows_gwasA/2)), rep(2, ceiling(n_rows_gwasA/2)))
gwasA$gwas.all2[2, "marker"] <- gwasA$gwas.all[1, "marker"]
####################################
test_that("manhattan.plot simple call", {
  expect_no_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = tag.table,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = tag.table,
      tag.repel = FALSE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all2,
      pvalue.thr = 5,
      tag.table = NULL,
      tag.repel = FALSE,
      point.colour = "red",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = "index",
      facet.main = "index",
      facet.secondary = NULL,
      facet.dim = c(20, 10),
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all2,
      pvalue.thr = 0.0005,
      tag.table = NULL,
      tag.repel = FALSE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = "index",
      facet.main = "index",
      facet.secondary = NULL,
      facet.dim = c(2, 1),
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all2,
      pvalue.thr = 0.0005,
      tag.table = NULL,
      tag.repel = FALSE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = "index",
      facet.main = "index",
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_no_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = tag.table,
      tag.repel = TRUE,
      point.colour = "",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
})
####################################
test_that("test traps on manhattan.plot", {
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all[,2:9],
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_no_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all[1:10000, ],
      pvalue.thr = 0.0005,
      tag.table = tag.table,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = "",
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = "",
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = tag.table[, c("chrom", "pos"), drop = FALSE],
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = "",
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = TRUE,
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = FALSE,
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
    )
  )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = "",
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = "",
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = "",
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = "",
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = "",
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = TRUE,
      facet.secondary = NULL,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = TRUE,
      facet.dim = NULL,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
  expect_error(
    manhattan.plot(
      gwas.table = gwasA$gwas.all,
      pvalue.thr = 0.0005,
      tag.table = gwasA$gwas.sel,
      tag.repel = TRUE,
      point.colour = "chrom",
      point.size = 1,
      point.alpha = 0.80,
      collate = TRUE,
      padding = 3e6,
      gwas.index = NULL,
      facet.main = NULL,
      facet.secondary = NULL,
      facet.dim = TRUE,
      scales = "free_y",
      legend.position = "none",
      message = TRUE
      )
    )
})
