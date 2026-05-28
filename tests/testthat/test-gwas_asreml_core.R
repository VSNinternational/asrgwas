testthat::skip_if_not_installed("asreml")
suppressPackageStartupMessages(library(asreml))
asreml.options(trace = FALSE)
all_markers <<- colnames(geno.apricot)
set.seed(1001)
gwas.data <- NULL
sample_sizes <- c(1000, 2000, 5000)
for (n_markers in sample_sizes) {
  sample_markers <- sample(all_markers, n_markers)
  new_marker <- all_markers[all_markers %in% sample_markers]
  gwas.data_try <- tryCatch(
    suppressWarnings(pre.gwas(
      pheno.data = pheno.apricot,
      indiv = "Ind", resp = "Sucrose",
      geno.data = geno.apricot[, new_marker, drop = FALSE],
      map.data = map.apricot[map.apricot$marker %in% new_marker, ],
      maf = 0, heterozygosity = 0.9, Fis = 0.5,
      marker.callrate = 1, ind.callrate = 1, impute = TRUE)),
    error = function(e) NULL
  )
  if (!is.null(gwas.data_try)) {
    gwas.data <- gwas.data_try
    break
  }
}
if (is.null(gwas.data)) {
  stop("Unable to create a non-empty GWAS fixture from sampled markers.")
}
gwas.data <<- gwas.data
geno.data2 <- gwas.data$geno.data
geno.data2[1,1] <- NA
pheno.data2 <- gwas.data$pheno.data
pheno.data2$Sucrose[pheno.data2$Sucrose %in% sample(pheno.data2$Sucrose, 200)] <- NA
model <<- asreml(fixed = Sucrose~1+Ind,
                 residual=~idv(units),
                 na.action=list(x='include',y='include'),
                 data=gwas.data$pheno.data)
model_no_converg <- model
model_no_converg$converge <- FALSE
gwas.data$pheno.data$binom <- rbinom(nrow(gwas.data$pheno.data), 1, 0.5)
gwas.data$pheno.data$binom_error <- rbinom(nrow(gwas.data$pheno.data), 2, 0.5)
suppressWarnings(model2 <<- asreml(binom ~ Ind,
                  random = ~ Lots,
                  family = asr_binomial(),
                  data = gwas.data$pheno.data))
caller <<- "gwas"
gwas.data_Kinv_Ind2 <- gwas.data$Kinv
names(gwas.data_Kinv_Ind2) <- "Ind2"
####################################
test_that("gwas.asreml simple call", {
  expect_no_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = model,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
    )
  )
  expect_no_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "binom",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = NULL,
      randomf = "Lots",
      residual = "Lots",
      weights = NULL,
      family = "binomial",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = model2,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
    )
  )
  expect_warning(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = "Lots",
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = geno.data2,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
    )
  )
  expect_warning(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = "Lots",
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = geno.data2,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = 1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
    )
  )
  expect_no_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = NULL,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb")
  )
  expect_no_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = NULL,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb")
  )
  expect_no_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = NULL,
      map.data = NULL,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb")
  )
  expect_no_error(
    gwas.asreml(
      pheno.data = data.table::as.data.table(gwas.data$pheno.data),
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = NULL,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb")
  )
  expect_no_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = NULL,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb",
      P3D = FALSE
      )
  )
  expect_no_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data[-1,],
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = NULL,
      pvalue.thr = 0.0005,
      bonferroni = FALSE
      )
  )
})
####################################
test_that("test traps on gwas.asreml", {
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = NULL,
      randomf = NULL,
      residual = "Lots",
      mod = model_no_converg,
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      bonferroni = TRUE
      )
  )
  expect_error(
    gwas.asreml(
      pheno.data = pheno.data2,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = NULL,
      randomf = NULL,
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      P3D = TRUE
      )
  )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = c("Ind","Ind"),
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = NULL,
      map.data = NULL,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb"
      )
  )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      randomf = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = NULL,
      map.data = NULL,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb"
    )
  )
  expect_error(
    gwas.asreml(
      pheno.data = as.matrix(gwas.data$pheno.data),
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      family = "binomial",
      gen = "Ind",
      fixedf = NULL,
      randomf = NULL,
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      bonferroni = TRUE
      )
    )
  expect_error(
    ASRgwas::gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = NULL,
      randomf = NULL,
      residual = "Lots",
      Kinv = gwas.data_Kinv_Ind2,
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      bonferroni = TRUE,
      P3D = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = list(gwas.data$Kinv$Ind[-1,]),
      Q = gwas.data$Q,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb"
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = list(gwas.data$Kinv$Ind[1:3,]),
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = NULL,
      npc = 3,
      weights = "Lots",
      family = "gaussian",
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb")
  )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = NULL,
      npc = 3,
      weights = "Lots",
      family = "binomial",
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb")
  )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "binom_error",
      total = "Ind",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      family = "binomial"
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "binom_error",
      total = NULL,
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      family = "binomial"
    )
  )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      fixedf = "Lots",
      residual = "Lots",
      Kinv = gwas.data$Kinv,
      Q = NULL,
      npc = 3,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data[,2:3],
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      workspace = "2Gb")
  )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = -1,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = -0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = c("Lots", "Years"),
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = c("Lots", "Years"),
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = c("Sucrose", "Glucose"),
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_error(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = -3,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
  expect_warning(
    gwas.asreml(
      pheno.data = gwas.data$pheno.data,
      resp = "Sucrose",
      gen = "Ind",
      Kinv = gwas.data$Kinv,
      Q = gwas.data$Q,
      npc = 1000,
      cov = NULL,
      fixedf = "Lots",
      randomf = NULL,
      residual = "Lots",
      weights = NULL,
      family = "gaussian",
      dispersion = 1,
      total = NULL,
      workspace = "1Gb",
      mod = NULL,
      geno.data = gwas.data$geno.data,
      map.data = gwas.data$map.data,
      pvalue.thr = 0.0005,
      bonferroni = FALSE,
      P3D = TRUE,
      maxiter.update = 2L,
      inverse.update = -1L,
      threads = detectCores() - 1,
      obj.parallel = NULL,
      message = TRUE
      )
    )
})
