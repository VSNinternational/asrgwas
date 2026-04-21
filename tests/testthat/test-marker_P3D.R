brute_force_p3d <- function(y, X, geno.data, V, nedf) {
  count_NA.marker <- colSums(is.na(geno.data))
  effect <- effectvar <- numeric(ncol(geno.data))
  for (i in seq_len(ncol(geno.data))) {
    keep <- is.finite(geno.data[, i])
    marker <- geno.data[keep, i]
    cur_X <- cbind(marker, X[keep, , drop = FALSE])
    cur_Y <- matrix(y[keep], ncol = 1)
    cur_Vinv <- Matrix::chol2inv(Matrix::chol(V[keep, keep, drop = FALSE]))
    XVX <- crossprod(cur_X, cur_Vinv %*% cur_X)
    XVy <- crossprod(cur_X, cur_Vinv %*% cur_Y)
    XVX_i <- solve(XVX)
    cur_beta <- XVX_i %*% XVy
    effect[i] <- cur_beta[1, 1]
    effectvar[i] <- XVX_i[1, 1]
  }
  std.error <- sqrt(effectvar)
  z.ratio <- effect / std.error
  p.value <- exp(
    log(2) +
      pt(abs(z.ratio), lower.tail = FALSE, log.p = TRUE, df = (nedf - count_NA.marker - 1))
  )
  list(
    effect = effect,
    std.error = std.error,
    z.ratio = z.ratio,
    p.value = p.value
  )
}
test_that("P3D exact mixed-missing path matches brute-force GLS", {
  y <- c(1.3, -0.2, 0.7, 1.9, -1.1, 0.4)
  X <- cbind(
    intercept = 1,
    covariate = c(-2, -1, 0, 1, 2, 3)
  )
  geno.data <- cbind(
    marker_1 = c(0, 1, 2, 0, 1, 2),
    marker_2 = c(1, NA, 0, 1, 2, 0),
    marker_3 = c(2, 1, 0, NA, 1, 2),
    marker_4 = c(0, 0, 1, 1, 2, NA)
  )
  V <- toeplitz(c(1, 0.35, 0.18, 0.09, 0.04, 0.02))
  nedf <- 24
  expected <- brute_force_p3d(
    y = y,
    X = X,
    geno.data = geno.data,
    V = V,
    nedf = nedf
  )
  observed <- ASRgwas:::P3D(
    y = y,
    X = X,
    geno.data = geno.data,
    V = V,
    nedf = nedf,
    inverse.update = -1,
    threads = 1,
    message = FALSE
  )
  expect_equal(observed$effect, expected$effect, tolerance = 1e-10)
  expect_equal(observed$std.error, expected$std.error, tolerance = 1e-10)
  expect_equal(observed$z.ratio, expected$z.ratio, tolerance = 1e-10)
  expect_equal(observed$p.value, expected$p.value, tolerance = 1e-10)
})
