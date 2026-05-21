#' Implements Population Parameters Previously Determined (P3D) algorithm
#'
#' @param y A vector with complete (no missing) centred phenotypic values of length \eqn{n} (default = \code{NULL}).
#' @param X A design matrix of fixed effects of size \eqn{n \times f}, where \eqn{n} is the
#' number of samples (as in \code{y}) and \eqn{f} is the number of fixed effects to be fitted along
#' with markers (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p},
#' with \eqn{n} individuals and \eqn{p} markers. Individual and marker names are
#' assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as: 0, 1, 2 (default = \code{NULL}).
#' @param V A matrix of size \eqn{n \times n} with the variance-covariance of \eqn{y} (default = \code{NULL}).
#' @param nedf The number of degrees of freedom from the base model. This is obtained from
#' \code{mod$nedf} (default = \code{NULL}).
#' @param inverse.update An integer indicating the number of Schulz-type iterative updates to be applied to
#' the inverse of \eqn{\boldsymbol{V}} = \emph{var(y)} if \emph{geno.data}
#' has missing values when running a Gaussian model under P3D.
#' The larger the value the better the approximation (under certain conditions).
#' If \code{inverse.update = -1} then the Woodbury's method is used instead.
#' If \code{inverse.update = 0} is used, then Schulz's procedure is performed with no iterations (no update),
#' been less precise but faster to run (default = \code{-1}).
#' @param threads An integer with the number of threads to be used in parallel processing
#' (default = \code{1}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with \code{effect}, \code{z.ratio} and \code{p.value} of each marker present
#' in \emph{geno.data}.
#'
#' @keywords internal
P3D <- function(y = NULL, X = NULL, geno.data = NULL, V = NULL, nedf = NULL,
                inverse.update = 0, threads = 1, message = TRUE){
  count.NA.marker <- colSums(is.na(geno.data))
  any.miss.marker <- any(count.NA.marker > 0)
  if (message){
    message("Inverting var(y) using matrix equations.")
  }
  try(Vinv <- Matrix::chol2inv(Matrix::chol(V)))
  if (!exists("Vinv")){
    stop("The inverse of var(y) cannot be obtained for the provided base model.")
  }
  y.mat <- as.matrix(y)
  Vinv.mat <- as.matrix(Vinv)
  if (message){
    message("Iterating through markers...")
  }
  time.start <- proc.time()["elapsed"]
  if (!any.miss.marker){
    effect.n.var <-
      P3D_cpp(
        Y = y.mat, M = geno.data, X = X,
        Vinv = Vinv.mat, ncores = threads)
    effect <- effect.n.var[, 1]
    effectvar <- effect.n.var[, 2]
  }
  if (any.miss.marker & inverse.update == -1){
    effect <- effectvar <- rep(NA_real_, ncol(geno.data))
    complete.marker <- count.NA.marker == 0
    missing.marker <- !complete.marker
    if (any(complete.marker)) {
      effect.n.var <-
        P3D_cpp(
          Y = y.mat,
          M = geno.data[, complete.marker, drop = FALSE],
          X = X,
          Vinv = Vinv.mat,
          ncores = threads
        )
      effect[complete.marker] <- effect.n.var[, 1]
      effectvar[complete.marker] <- effect.n.var[, 2]
    }
    if (any(missing.marker)) {
      effect.n.var <-
        P3D_woodbury_cpp(
          Y = y.mat,
          M = geno.data[, missing.marker, drop = FALSE],
          X = X,
          Vinv = Vinv.mat,
          ncores = threads
        )
      effect[missing.marker] <- effect.n.var[, 1]
      effectvar[missing.marker] <- effect.n.var[, 2]
    }
  }
  if (any.miss.marker & inverse.update != -1){
    y <- do.call("cbind", replicate(ncol(geno.data), y, simplify = FALSE))
    iter <- 1:ncol(geno.data)
    for (i in iter) {y[is.na(geno.data[, i]), i] <- 0}
    X <- scale(X, center = TRUE, scale = FALSE)
    geno.data <- scale(geno.data, center = TRUE, scale = FALSE)
    effect.n.var <-
      P3D_schulz_cpp(
        Y = as.matrix(y), M = geno.data, X = X,
        V = as.matrix(V), Vinv = Vinv.mat,
        niter = inverse.update, ncores = threads)
    effect <- effect.n.var[, 1]
    effectvar <- effect.n.var[, 2]
  }
  time.finish <- proc.time()["elapsed"] - time.start
  std.error <- sqrt(effectvar)
  z.ratio <- effect/std.error
  p.value <- exp(
    log(2) +
      pt(abs(z.ratio), lower.tail = FALSE, log.p = TRUE, df = (nedf - count.NA.marker - 1))
  )
  return(list(effect = effect, std.error = std.error,
              z.ratio = z.ratio, p.value = p.value,
              time.taken = time.finish))
}
