#' Selects final set of GWAS markers using backward elimination
#'
#' @description
#' The backward elimination technique starts by adding all markers to the base model as covariates
#' and obtaining their \eqn{p-value}. The marker with the highest \eqn{p-value} is eliminated and
#' the model is re-fit. The process is repeated until all retained markers have a \eqn{p-value}
#' smaller than the provided threshold.
#'
#' @param gwas.object A (list) with all objects obtained
#' after running the function \link{gwas.asreml} (default = \code{NULL}).
#' @param geno.data.sel A matrix with (a subset of) marker data of form \eqn{n \times s}, with \eqn{n} individuals and
#' \eqn{s} "significant" markers to be evaluated in the backward selection process. Individuals and marker names
#' are assigned to \code{rownames} and \code{colnames}, respectively.
#' Markers can have any coding (\emph{e.g.}, additive or dominant) as they will be fitted as covariates in the model (default = \code{NULL}).
#' @param ref.vc A numeric (vector) containing the indices (or positions) of the variance components
#' (from \code{gwas.object$mod$vparameters}) to be used as the reference for the variance explained
#' by the markers (usually the genetic variance) (default = \code{1}).
#' @param pvalue.thr A numeric value with the corresponding \eqn{p-value}
#' threshold to identify significant markers for the elimination process.
#' Markers will be eliminated until all have a \eqn{p-value} smaller
#' than the value provided (default = \code{0.1})
#' @param maxiter.update A numeric value indicating the number of iterations to be
#' used in \link[asreml]{update.asreml} (default = \code{3L}).
#' @param workspace Specifies the workspace to be used by the \code{asreml} function
#' (default = \code{"1Gb"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return
#' A list containing:
#' \itemize{
#'   \item \code{call}: An object of class \code{call} with the code used to fit \code{mod}.
#'   \item \code{gwas.sel}: A vector with the names of the final marker set with \code{p-value < pvalue.thr}.
#'   \item \code{mod}: The GWAS model fitted with the retained markers. Contains the same set of columns as \code{gwas.sel} from \link{gwas.asreml}.
#'   \item \code{aov}: The ANOVA table of the updated model with the final marker set.
#'   \item \code{expl.var}: The percentage of the genetic variance explained by the final marker set.
#'   \item \code{collinear.markers} [returned if collinearity is an issue]: A data frame with the names and correlation values of highly collinear/correlated markers that might affect model fitting.
#' }
#'
#' @details
#' If there is a group of markers with high collinearity (\emph{i.e.}, correlation >= 0.95)
#' then the process stops and the list of markers that are highly correlated
#' is reported. Here, one or more of these markers needs to be eliminated
#' for the backward elimination to proceed.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare Apricot dataset.
#' gwas.data <- pre.gwas(
#'   pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
#'   geno.data = geno.apricot, map.data = map.apricot, maf = 0.05,
#'   marker.callrate = 0.40, impute = TRUE)
#'
#' # Run GWAS.
#' gwasA <- gwas.asreml(
#'   pheno.data = gwas.data$pheno.data, resp = "Sucrose", gen = "Ind",
#'   fixedf = "Lots", residual = "Lots",
#'   Kinv = gwas.data$Kinv, Q = gwas.data$Q, npc = 3,
#'   geno.data = gwas.data$geno.data, map.data = gwas.data$map.data,
#'   pvalue.thr = 0.0005, bonferroni = FALSE, workspace = "2Gb")
#'
#' # Get name of significant markers.
#' sign.markers <- gwasA$gwas.sel$marker
#' geno.data.sel <- gwas.data$geno.data[, sign.markers]
#'
#' # Identify the index of genetic variance components.
#' gwasA$mod$vparameters
#'
#' # Run of marker selection (this will return only the correlations).
#' set <- select.marker(
#'   gwas.object = gwasA, geno.data.sel = geno.data.sel,
#'   ref.vc = 1, pvalue.thr = 0.1)
#' set$call
#' set$aov
#' }
select.marker <- function(
    gwas.object = NULL,
    geno.data.sel = NULL,
    ref.vc = 1,
    pvalue.thr = .1,
    maxiter.update = 5,
    workspace = "1Gb",
    message = TRUE
){
  .argument_class(.data = gwas.object, .class = "gwas.asreml")
  .argument_class(.data = geno.data.sel, .class = c("matrix", "array"))
  if(is.null(rownames(geno.data.sel))){
    stop("Individual names not assigned to rows of `geno.data,sel`.")
  }
  if(is.null(colnames(geno.data.sel))){
    stop("Marker names not assigned to columns of `geno.data.sel`.")
  }
  if (pvalue.thr <= 0 || pvalue.thr >= 1){
    stop("The `pvalue.thr` argument should follow the condition 0 < pvalue.thr < 1.")
  }
  if (maxiter.update < 0){
    stop("The `maxiter.update` aregument should follow the condition 0 <= maxiter.update.")
  }
  if (!all(ref.vc %in% 1:length(gwas.object$mod$vparameters))){
    stop("The provided `ref.vc` does not match indices of `gwas.object$mod$vparameters`.",
         " Maximum index is ", length(gwas.object$mod$vparameters), ".")
  }
  .argument_class(.data = message, .class = "logical")
  if (message){
    message(col_blue('\nChecking datasets.'))
  }
  cor.report <- correlation.report_(matrix_ = geno.data.sel, cor.thr_ = .95)
  if (!is.null(cor.report)) {
    return(list(collinear.markers = cor.report))
  }
  mod <- gwas.object$mod
  pheno.data <- mod$mf
  resp <- attributes(pheno.data)$traits$lhs
  gen.termvar <- names(mod$vparameters[ref.vc])
  family <- attributes(pheno.data)$traits$family
  if (message){
    message("Assuming `", names(pheno.data)[1], "` (first column of `pheno.data`) as genotype identifier.")
  }
  indiv <- names(pheno.data)[1]
  if (!identical(as.character(pheno.data[[indiv]]), rownames(geno.data.sel))){
    geno.data.sel <- geno.data.sel[as.character(pheno.data[[indiv]]), ]
  }
  try.variables <- colnames(geno.data.sel)
  data = cbind.data.frame(pheno.data, geno.data.sel)
  rm(geno.data.sel)
  if (message){
    message(col_blue('\nPerforming backward elimination of markers.'))
    message("This function calls `asreml` several times which may take a while.")
  }
  bw.elim <- backward.elimination(
    mod = mod,
    data = data,
    try.variables = try.variables,
    pvalue.thr = pvalue.thr,
    maxiter.update = maxiter.update,
    workspace = workspace,
    message = message)
  if (length(bw.elim$kept.variables) == 0){
    if (message){
      message("No markers kept in GWAS model according to p-value threshold of ", pvalue.thr, ".")
    }
    bw.elim$mod <- gwas.object$mod
    bw.elim$aov <- gwas.object$aov
    gwas.sel <- NULL
    expl.var <- NULL
  } else {
    gwas.sel <- gwas.object$gwas.all[gwas.object$gwas.all$marker %in% bw.elim$kept.variables,
                                     c("marker", "chrom", "pos", "maf")]
    gwas.sel <- cbind(
      gwas.sel,
      summary(bw.elim$mod, coef = TRUE)$coef.fixed[gwas.sel$marker,])
    names(gwas.sel)[5] <- "effect"
    gwas.sel$p.value <- bw.elim$aov[gwas.sel$marker, "Pr"]
    if (family == "gaussian") {
      gwas.sel$expl.var <- 100 * (2 * gwas.sel$maf * (1 - gwas.sel$maf) * (gwas.sel$effect^2)) / var(data[[resp]], na.rm = TRUE)
    }
    rownames(gwas.sel) <- NULL
    if (message){
      message(col_blue(paste0('\nSummary report.')))
      message("A total of ", length(bw.elim$kept.variables),
              " marker(s) kept according to threshold of ", pvalue.thr, ".")
    }
    if (!is.null(ref.vc) & family == "gaussian"){
      expl.var <- round((1 - bw.elim$mod$vparameters[gen.termvar] /
                           mod$vparameters[gen.termvar]) * 100, 4)
      expl.var[expl.var < 0] <- 0
      if (message){
        message(paste0("The set of kept markers explain ", expl.var,
                       "% of the variance of ", gen.termvar, ".", collapse = "\n"))
      }
      expl.var <- as.data.frame(expl.var)
      names(expl.var) <- "expl.var"
    } else {
      expl.var <- NULL
    }
  }
  bw.elim$mod$call <- fix.call.param_(mod_ = bw.elim$mod)
  return(
    list(
      call = bw.elim$mod$call,
      gwas.sel = gwas.sel,
      mod = bw.elim$mod,
      aov = bw.elim$aov,
      expl.var = expl.var))
}
