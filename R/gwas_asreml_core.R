#' Directs data frames and matrices to the respective GWAS method
#'
#' @description
#' Passes the data frames and matrices to the available GWAS model fit functions.
#' Currently methods: \code{P3D} and \code{marker.update} are available.
#'
#' @param mod A pre-fit \pkg{asreml} model object. Must follow the conditions of the function
#' (default = \code{NULL}).
#' @param weights A character with name of the numeric variable in \code{pheno.data} with the
#' weights of each prediction (\emph{i.e.}, mean estimate) (default = \code{NULL}).
#' @param map.data A data frame with each marker's name, chromosome, and position.
#' Variable names \strong{must} be "marker", "chrom", and "pos" (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p},
#' with \eqn{n} individuals and \eqn{p} markers. Individual and marker names are
#' assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as: 0, 1, 2 (default = \code{NULL}).
#' @param Kinv A list of matrices representing the \strong{inverse} of the genomic relationship
#' matrix \eqn{\boldsymbol{Kinv}}. The matrices must be in sparse form (default = \code{NULL}).
#' @param pvalue.thr A numeric value with \eqn{p-value} threshold to identify significant markers (default = \code{5e-6})
#' @param P3D If \code{TRUE} the "Population Parameters Previously Determined" algorithm will be used.
#' If \code{FALSE} the variance components are allowed to vary while fitting each marker individually.
#' In setting \code{P3D = TRUE} there is a considerable speed gain for \code{family = "gaussian"}, but not for
#' \code{family = "binomial"} (default = \code{TRUE}).
#' @param user.mod Set to \code{TRUE} if the base model has been calculated outside \link{gwas.asreml}
#' (default = \code{NULL}).
#' @param maxiter.update A numeric value indicating the number of iterations to be
#' used in \link[asreml]{update.asreml} under the non-P3D method (default = \code{1}).
#' @param threads An integer with the number of threads to be used in parallel processing
#' (default = \code{1}).
#' @param obj.parallel A character vector indicating which objects passed to the function that must
#' be exported to clusters (\emph{e.g.}, objects passed to \code{Kinv} argument). If nothing is passed, the
#' algorithm will try to identify the required objects (default = \code{NULL}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with several outputs from the GWAS model fit that will be passed to \link{gwas.asreml}:
#'
#' \itemize{
#'   \item \code{call}: The \pkg{asreml} model call.
#'   \item \code{gwas.all}: A data frame with map and marker statistics of all tested markers.
#'   \item \code{gwas.sel}: A data frame with map and marker statistics of significant markers.
#'   \item \code{mod}: The \pkg{asreml} object GWAS model fitted with the provided data.
#' }
#'
#' @keywords internal
gwas.asreml.core <- function(
    mod = NULL,
    weights = NULL,
    map.data = NULL,
    geno.data = NULL,
    Kinv = NULL,
    pvalue.thr = 5e-6,
    P3D = TRUE,
    user.mod = NULL,
    maxiter.update = 1,
    cpp = NULL,
    inverse.update = NULL,
    threads = 1,
    obj.parallel = NULL,
    message = TRUE
) {
  data <- mod$mf
  resp <- attributes(data)$traits$lhs
  family <- attributes(mod$mf)$traits$family
  if (message){
    message(col_blue("\nStarting Genome-Wide Association model fit."))
    message("A total of ", ncol(geno.data), " markers will be evaluated on ", nrow(geno.data), " samples.")
  }
  if (P3D & family == "gaussian" & !user.mod) {
    if (message){
      message("Initiating P3D estimation of marker effects.")
      message("Estimating V = var(y) from base model fit.")
    }
    K <- lapply(Kinv, function(k) ASRgenomics::sparse2full(k))
    K <- lapply(K, function(k) chol2inv(chol(k)))
    K <- lapply(K, function(k){
      rownames(k) <- colnames(k) <- attr(Kinv[[1]], "rowNames")
      return(k)
    }
    )
    names(K) <- names(Kinv)
    V <- v.asreml(asreml.obj = mod, Ks = K, weights = weights)
    rm(K)
    n.fix.eff <- nrow(mod$coefficients$fixed)
    if (message){
      message("Constructing design matrices (with a total of ", n.fix.eff, " fixed effect levels identified).")
    }
    X <- as.matrix(mod$design)[, seq(length.out = n.fix.eff), drop = FALSE]
    y <- data[[resp]]
    gwas.all <- P3D(y = y, X = X, geno.data = geno.data, V = V,
                    nedf = mod$nedf, inverse.update = inverse.update,
                    threads = threads, message = message)
  }
  else {
    if (message){
      message("Initiating estimation of marker effects using asreml `update` method.")
      if (threads > 1){
          message("Note: maximum `threads` is constrained to the number of concurrent sessions allowed by the ASReml-R license.")
      }
      geno.data.time.sample <- 1:(threads * 5)
      if (max(geno.data.time.sample) < ncol(geno.data)){
        runtime <-
          mc.marker.update(
            mod = mod, geno.data = geno.data[, geno.data.time.sample],
            maxiter.update = maxiter.update,
            P3D = P3D, threads = threads,
            obj.parallel = obj.parallel,
            message = FALSE)$time.taken
        message("Expected runtime: ",
                (ncol(geno.data) * as.numeric(runtime) / max(geno.data.time.sample) / 60) |> round(2), " min.")
      }
    }
    gwas.all <-
      mc.marker.update(
        mod = mod, geno.data = geno.data,
        maxiter.update = maxiter.update,
        P3D = P3D, threads = threads,
        obj.parallel = obj.parallel,
        message = message)
  }
  if (message){
    message("Actual runtime: ",
            (as.numeric(gwas.all$time.taken) / 60) |> round(2), " min.")
  }
  if (family == "gaussian") {
    expl.var <- 100 * (2 * map.data$maf * (1 - map.data$maf) * (gwas.all$effect^2)) / var(data[[resp]], na.rm = TRUE)
  }
  gwas.all <- data.frame(
    marker = map.data$marker,
    chrom = map.data$chrom,
    pos = map.data$pos,
    maf = map.data$maf,
    effect = gwas.all$effect,
    std.error = gwas.all$std.error,
    z.ratio = gwas.all$z.ratio,
    p.value = gwas.all$p.value
  )
  rownames(gwas.all) <- NULL
  if (family == "gaussian") {
    gwas.all <- cbind.data.frame(gwas.all, expl.var)
  }
  gwas.sel <- gwas.all[
    !is.na(gwas.all$p.value) &
      gwas.all$p.value <= pvalue.thr,]
  if (message){
    message(
      col_blue("\nA total of ", nrow(gwas.sel),
               " markers identified based on a p-value threshold of ",
               pvalue.thr, "."))
  }
  if (nrow(gwas.sel) == 0){
    gwas.sel = NULL
  }
  if (message & !is.null(gwas.sel)){
    fdr <- nrow(gwas.all) * pvalue.thr / nrow(gwas.sel)
    message("The observed FDR for the set of significant markers is: ", (fdr * 100) |> round(4), "%.\n")
  }
  return(list(
    call = mod$call,
    gwas.all = gwas.all,
    gwas.sel = gwas.sel,
    mod = mod
  ))
}
