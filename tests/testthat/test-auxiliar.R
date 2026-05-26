G <- G.matrix(
  M = ASRgenomics::geno.apple,
  method = "VanRaden",
  sparseform = FALSE
)$G
new_G <- G[1:5, 1:5]
new_G_init <- new_G
new_G_init[1, 1] <- NA
woodbury_base <- matrix(
  c(4, 1, 2, 0, 1,
    1, 5, 1, 1, 0,
    2, 1, 6, 1, 1,
    0, 1, 1, 4, 1,
    1, 0, 1, 1, 3),
  nrow = 5, byrow = TRUE
)
woodbury_inv <- solve(woodbury_base)
####################################
test_that("moore.penrose simple call", {
  expect_no_error(
    moore.penrose(
      G = new_G,
      eig.tol = 1e-6
    )
  )
  expect_no_error(
    moore.penrose(
      G = new_G,
      eig.tol = 2
    )
  )
})
####################################
test_that("woodbury simple call", {
  expect_no_error(
    woodbury(
      X = new_G,
      na = c(1, 2)
    )
  )
})
####################################
test_that("woodbury_cpp matches exact row-column elimination", {
  dropped <- c(1L, 3L)
  kept <- setdiff(seq_len(nrow(woodbury_base)), dropped)
  expect_equal(
    ASRgwas:::woodbury_cpp(
      X = woodbury_inv,
      na = dropped - 1L,
      nona = kept - 1L
    ),
    solve(woodbury_base[kept, kept, drop = FALSE]),
    tolerance = 1e-10
  )
  expect_equal(
    ASRgwas:::woodbury_cpp(
      X = woodbury_inv,
      na = dropped - 1L,
      nona = kept - 1L
    ),
    woodbury(
      X = woodbury_inv,
      na = dropped
    ),
    tolerance = 1e-10
  )
})
####################################
test_that("schulz simple call", {
  expect_no_error(
    schulz(
      X = new_G,
      Xinv.init = new_G_init,
      na = c(1, 2),
      niter = 2
    )
  )
})
