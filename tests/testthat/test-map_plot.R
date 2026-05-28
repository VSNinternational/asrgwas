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
####################################
test_that("map.plot simple call", {
  expect_no_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker)
  )
  expect_no_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
})
####################################
test_that("test traps on map.plot", {
  expect_error(
    map.plot(map.data = gwasA$gwas.all[ , 2:9],
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = 1,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
    )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = "",
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = 1,
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = "1.8",
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = 1,
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = "0.05",
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = "0.08",
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = 1,
             chrom.width = 4,
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = "4",
             chrom.contour = "#DEDEDE",
             message = TRUE)
  )
  expect_error(
    map.plot(map.data = gwasA$gwas.all,
             tag.markers = gwasA$gwas.sel$marker,
             tag.repel = FALSE,
             tag.colour = "black",
             tag.size = 1.8,
             marker.colour = "#105E26",
             marker.width = 0.05,
             marker.alpha = 0.08,
             chrom.colour = "#DEDEDE",
             chrom.width = 4,
             chrom.contour = 1,
             message = TRUE)
  )
})
