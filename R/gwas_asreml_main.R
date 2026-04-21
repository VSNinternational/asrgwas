#' Fits Genome-wide Association Studies (GWAS) models using asreml
#'
#' @description
#' Main function of \pkg{ASRgwas} that extends
#' the capabilities of \pkg{asreml} to perform Genome-wide Association Studies.
#'
#' @param pheno.data A data frame with all relevant columns (fphenotypic response,
#' factors and covariates) to be used (default = \code{NULL}).
#' @param resp A character string of the name of the numeric variable in \code{pheno.data}
#' that identifies the response variable (default = \code{NULL}).
#' @param gen A character (vector) with the names of the columns in \code{pheno.data}
#' that identifies the genotypes for the "genetic" part of the model (default = \code{NULL}).
#' @param Kinv A list of matrices representing the \strong{inverse} of the genomic relationship
#' matrix \eqn{\boldsymbol{Kinv}}. The matrices must be in sparse form (default = \code{NULL}).
#' @param Q A matrix of size \eqn{n} x \eqn{nv}, with \eqn{nv} population structure-related covariates (\emph{e.g.},
#' principal component scores). The number of components to be included in the fitted model
#' is determined by \code{npc} (default = \code{NULL}).
#' @param npc The number of components (columns) to be used from the \eqn{\boldsymbol{Q}} matrix
#' (default = \code{10L}).
#' @param cov A character vector with the names of the numeric variables in \code{pheno.data}
#' to be considered as covariates in the GWAS model (default = \code{NULL}).
#' @param fixedf A character vector with the names of the factors in \code{pheno.data}
#' to be considered as fixed effects in the GWAS model. Interactions can be added as \code{"factor1:factor2"}
#' (default = \code{NULL}).
#' @param randomf A character vector with the names of the factors in \code{pheno.data}
#' to be considered as random effects in the GWAS model. Interactions can be added as
#' \code{"factor1:factor2"}. Interactions with \code{gen} factor can also
#' be added here (\emph{e.g.}, \code{"gen:factor1"}) (random by interaction with \code{gen}).
#' Genotypes should be added to the argument \code{gen}, not the \code{randomf} (default = \code{NULL}).
#' @param weights A character stirng of the name of the numeric variable in \code{pheno.data} that
#' specifics the weights of each prediction (\emph{i.e.}, mean estimate) (default = \code{NULL}).
#' @param residual A character string of the name of the variable in \code{pheno.data} to be used as a
#' conditional factor for specifying a heterogeneous residual structure according to the levels of this factor.
#' The \code{pheno.data} will be internally ordered by this factor (default = \code{NULL}).
#' @param family A character indicating the family of the distribution for \code{resp}.
#' Options are: \code{"gaussian"} and \code{"binomial"} (default = \code{"gaussian"}).
#' @param dispersion A numeric value of the dispersion parameter
#' (used if \code{family = "binomial"}) (default = \code{1}).
#' @param total A character string of the name of the numeric variable in \code{pheno.data} that
#' indicates the total counts (if \code{family = "binomial"}). If \code{NULL}, count is taken as 1
#' (default = \code{NULL}).
#' @param workspace Specifies the workspace to be used by the \code{asreml} function.
#' Note that large workspaces will result in much longer processing times if \code{P3D} is false
#' (default = \code{"1Gb"}).
#' @param mod A pre-fitted \pkg{asreml} model object. Providing a
#' pre-fitted model works for Gaussian
#' data if \code{P3D = FALSE} and binomial data (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p},
#' with \eqn{n} individuals and \eqn{p} markers. Individual and marker names are
#' assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as: 0, 1, 2 (default = \code{NULL}).
#' @param map.data A data frame with the marker names, chromosomes, and positions.
#' Variable names \strong{must} be: "marker", "chrom", and "pos" (default = \code{NULL}).
#' @param pvalue.thr A numeric value with the \eqn{p-value} threshold to identify significant markers (default = \code{0.05})
#' @param bonferroni If \code{TRUE} the Bonferroni correction is used with expression \eqn{1 - (1 - pvalue.thr)^{1/p} \sim pvalue.thr/p},
#' where \eqn{pvalue.thr} is based on the argument (\code{pvalue.thr}) and \eqn{p} is the number of markers (default = \code{TRUE}).
#' @param P3D If \code{TRUE} the "Population Parameters Previously Determined" algorithm will be used.
#' If \code{FALSE} the variance components are allowed to vary while fitting each marker individually.
#' Setting \code{P3D = TRUE} results in considerable speed gain for \code{family = "gaussian"}, but not for
#' \code{family = "binomial"} (default = \code{TRUE}).
#' @param maxiter.update A numeric value indicating the number of iterations to be
#' used in \code{update.asreml} method if requested (default = \code{2}).
#' @param inverse.update An integer indicating the number of Schulz-type iterative updates to be applied to
#' the inverse of \eqn{\boldsymbol{V}} (phenotypic variance matrix) if \code{geno.data}
#' has missing values when running a Gaussian model under P3D.
#' The larger the value the better the approximation (under certain conditions).
#' If \code{inverse.update = -1} then the Woodbury's method is used instead.
#' If \code{inverse.update = 0} is used, then Schulz's procedure is performed with no iterations (no update),
#' which is less precise but faster to run. See details for more information (default = \code{-1}).
#' @param threads An integer with the number of thread to be used in parallel processing
#' (default = \code{parallel::detectCores() - 1}; \emph{i.e.}, total number of threads available minus one).
#' @param obj.parallel A character vector indicating which objects to be passed to the function that must
#' be exported to clusters (\emph{e.g.} objects passed to \code{Kinv} argument). If nothing is passed, the
#' algorithm will try to identify the required objects (default = \code{NULL}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details
#'
#' Some relevant considerations are presented below.
#'
#' Marker only models (no \eqn{\boldsymbol{Q}} or \eqn{\boldsymbol{K}}) can be fitted by specifying
#' \code{gen}, \code{Q}, and \code{Kinv} as \code{NULL}.
#'
#' The factors passed to \code{gen} are considered random, as are the ones passed via \code{randomf}.
#' The arguments are separated so it is clear which ones should be on the numerator (\code{gen})
#' and denominator (\code{gen + randomf}) of the heritability (or repeatability) formula.
#' Interactions with \code{gen} factors should be included in \code{randomf} and will be
#' always added to the denominator. No interactions between \code{gen} factors are allowed.
#'
#' If there is more than one type of genetic effect to be modeled (but associated
#' with different genomic matrices),
#' the column with genotype identifications in \code{pheno.data} should be copied with a
#' different name so there is one for each genetic effect considered.
#'
#' If relationship matrices are associated to \code{gen} factors,
#' they should be passed using \code{Kinv}. The name of objects in \code{Kinv} \strong{must}
#' match the values passed to \code{gen}. For example, to model additive and dominance effects,
#' the following code can be used: \code{gen = c("genA", "genD")}, where \code{"genA"} and
#' \code{"genD"} are columns in \code{pheno.data} with genotype identifiers;
#' and \code{Kinv = list("genA" = KinvA, "genD" = KinvD)}, where \code{KinvA} and
#' \code{KinvD} are inverse matrices of additive and dominance relationships in sparse form.
#' If a factor is passed to \code{gen} but not to \code{Kinv},
#' it is regarded as having independent levels (\emph{i.e.}, \code{id/idv} variance structure).
#' Also, all factors names passed to the command \code{Kinv} \strong{must} be present in
#' the command \code{gen}.
#'
#' Row and column names are checked for matching between and within lists, and between datasets.
#' In the case of the sparse form \code{Kinv}, attribute names \code{rowNames} and \code{colNames}
#' are used to perform this verification.
#'
#' Any number of fixed (factors and covariates) and random effects can be passed,
#' but only single trait analysis has been implemented so far.
#'
#' Missing values in \code{resp} and \code{weights} are allowed. They are dealt with by
#' eliminating the observations (rows) prior to the function call.
#'
#' Missing values on the marker are allowed in both \code{P3D = TRUE} and \code{P3D = FALSE}.
#' If \code{P3D = FALSE}, \pkg{asreml} deals with the missing data.
#' Two algorithms (Schulz and Woodbury) have been implemented to deal
#' with missing data by adjusting the inverse of the \eqn{\boldsymbol{V}}
#' (phenotypic variance) matrix if \code{P3D = TRUE}.
#' It is likely that Woodbury's method will be faster and more precise than
#' Schulz's method with one or more iterations. In addition, Schulz's method
#' is highly dependent on matrix stability, while Woodbury's method is more robust.
#' Note that the default is Woodbury's method.
#'
#' The variance explained by each marker present in \code{gwas.all/gwas.sel}
#' is calculated via \eqn{expl.var = 100 \times [2 \times maf \times (1-maf) \times effect^2] / var(resp)}.
#' The MAF (\code{maf}; and consequently the variance explained \code{expl.var}) are only correct
#' for \code{geno.data} with additive genotype coding.
#'
#' The calculation of the \eqn{\boldsymbol{V}} (phenotypic variance) matrix required by Gaussian P3D
#' is based on the code previously implemented for the library \pkg{asremlPlus} (Brien, 2021).
#'
#' @return A list of class \code{gwas.asreml} with several outputs from a GWAS model fit:
#'
#' \itemize{
#'   \item \code{call}: An object of class \code{call} with the code used to fit \code{mod} in \pkg{asreml}.
#'   \item \code{gwas.all}: A data frame with map (\code{marker}, \code{chrom}, \code{pos}),
#'     minor allele frequency (\code{maf}), \code{effect}, \code{z.ratio}, \code{p.value},
#'     and explained variance (\code{expl.var}; see details) of all tested markers.
#'   \item \code{gwas.sel}: A subset of \code{gwas.all} with significant markers based on \code{pvalue.thr}
#'     and/or Bonferroni.
#'   \item \code{mod}: The \pkg{asreml} object base model fitted with the provided data.
#'   \item \code{aov}: The ANOVA table of the fitted base model.
#'   \item \code{heritability}: The estimated heritability value for \code{h2.vc} and
#'     \code{h2.pev} definitions.
#' }
#'
#' @export gwas.asreml
#'
#' @references
#' Brien C. 2021.
#' asremlPlus: Augments 'ASReml-R' in fitting mixed models and packages generally in exploring prediction differences.
#' R package version 4.3-31, <https://CRAN.R-project.org/package=asremlPlus>.
#'
#' Hager, W.W. 1989. Updating the inverse of a matrix.
#' \emph{SIAM Review} 31(2):221–239.
#'
#' Petković M.D. 2014.
#' Generalized Schultz iterative methods for the computation of outer inverses.
#' \emph{Computers & Mathematics with Applications} 67(10).
#'
#' Zhang, Z. et al. 2010. Mixed linear model approach adapted for genome-wide association studies.
#' \emph{Nature Genetics}, 42(4):355–360.
#'
#' @examples
#' \dontrun{
#' # Prepare Apricot dataset.
#' gwas.data <- pre.gwas(
#'   pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
#'   geno.data = geno.apricot, map.data = map.apricot,
#'   maf = 0.05, heterozygosity = 0.9, Fis = 0.5,
#'   marker.callrate = 0.40, impute = TRUE)
#'
#' # Fit model to see how good it is.
#' gwasA <- gwas.asreml(
#'   pheno.data = gwas.data$pheno.data,
#'   resp = "Sucrose",
#'   gen = "Ind",
#'   fixedf = "Lots",
#'   residual = "Lots",
#'   Kinv = gwas.data$Kinv,
#'   Q = gwas.data$Q,
#'   npc = 3)
#'
#' # Check the call.
#' gwasA$call
#'
#' # Check heritability/repeatability.
#' gwasA$heritability
#'
#' # Check the ANOVA table.
#' gwasA$aov
#'
#' # Fit the GWAS for that model.
#' gwasA <- gwas.asreml(
#'   pheno.data = gwas.data$pheno.data,
#'   resp = "Sucrose",
#'   gen = "Ind",
#'   fixedf = "Lots",
#'   residual = "Lots",
#'   Kinv = gwas.data$Kinv,
#'   Q = gwas.data$Q,
#'   npc = 3,
#'   geno.data = gwas.data$geno.data,
#'   map.data = gwas.data$map.data,
#'   pvalue.thr = 0.0005,
#'   bonferroni = FALSE,
#'   workspace = "2Gb")
#'
#' # The quantile-quantile plot.
#' qq.plot(gwas.table = gwasA$gwas.all)
#'
#' # The Manhattan plot.
#' manhattan.plot(gwas.table = gwasA$gwas.all, padding = 2e6)
#'
#' # Statistics for all markers.
#' head(gwasA$gwas.all)
#'
#' # Statistics for only significant markers.
#' head(gwasA$gwas.sel)
#' }
gwas.asreml <- function(
    pheno.data = NULL,
    resp = NULL,
    gen = NULL,
    Kinv = NULL,
    Q = NULL,
    npc = 10L,
    cov = NULL,
    fixedf = NULL,
    randomf = NULL,
    residual = NULL,
    weights = NULL,
    family = c('gaussian', 'binomial'),
    dispersion = 1,
    total = NULL,
    workspace = "1Gb",
    mod = NULL,
    geno.data = NULL,
    map.data = NULL,
    pvalue.thr = 0.05,
    bonferroni = TRUE,
    P3D = TRUE,
    maxiter.update = 2L,
    inverse.update = -1L,
    threads = detectCores() - 1,
    obj.parallel = NULL,
    message = TRUE
) {
  if (is.null(attributes(resp)$caller))
    caller <- "gwas"
  else if (attributes(resp)$caller == "gblup")
    caller <- "gblup"
  if (message){
    message(col_blue(paste0('\nChecking input.')))
  }
  if (is.null(mod)) {
    if (!is.null(Kinv)) {
      Kinv.name <- deparse(substitute(Kinv))
    }
  }
  if (!is.null(mod)) {
    user.mod <- TRUE
  } else{
    user.mod <- FALSE
  }
  if (user.mod){
    if (message){
      message("Model related arguments will be collected from `mod` and ignored if provided.")
    }
    .argument_class(.data = mod, .class = "asreml")
    if (!mod$converge){
      stop("The provided base model did not converge. ",
           "Please use `asreml::update.asreml()` to allow for more iterations.")
    }
    orig.data.name <- deparse(mod$call$data)
    input.from.asreml(mod = mod)
  }
  .argument_class(.data = message, .class = "logical")
  .argument_class(.data = P3D, .class = "logical")
  family <- match.arg(family)
  if (maxiter.update < 0){
    stop("The `maxiter.update` argument should follow the condition 0 <= maxiter.update.")
  }
  if (pvalue.thr < 0){
    stop("The `pvalue.thr` argument should follow the condition 0 <= pvalue.thr")
  }
  if (length(residual) > 1){
    stop("Only a single element should be passed to `residual`.")
  }
  if (length(weights) > 1){
    stop("Only a single element should be passed to `weights`.")
  }
  if (length(resp) > 1){
    stop("Only a single element should be passed to `resp`.")
  }
  if (!is.null(Q)){
    if (npc < 0){
      stop("The `npc` argument should be a positive number.")
    }
    if (npc > ncol(Q)){
      warning("The number of required components is larger than available. Hence, `npc` was set to `ncol(Q)`.")
      npc <- ncol(Q)
    }
    if (is.null(npc)){
      warning("As `Q` is provided, `npc` should be a positive number. Hence, `npc` was set to `ncol(Q)`.")
      npc <- ncol(Q)
    }
  }
  .argument_class(.data = pheno.data, .class = "data.frame")
  pheno.data <- droplevels(pheno.data)
  if (!is.null(geno.data)){
    .argument_class(.data = geno.data, .class = c("matrix", "array", "ddiMatrix"))
    if (is.null(map.data)){
      map.data <- ASRgenomics:::dummy.map_(marker.id = colnames(geno.data), message = message)
    } else {
      .argument_class(.data = map.data, .class = "data.frame")
    }
    if (!identical(names(map.data), c("marker", "chrom", "pos"))){
      stop("Names on `map.data` must be `marker', 'chrom', and 'pos'.")
    }
  }
  if (is.null(geno.data) & caller == "gwas"){
    if (message){
      message(col_yellow(paste0("\nObject `geno.data` NOT provided, fitting model only.")))
    }
  }
  if (!user.mod){
    if (!is.null(Kinv)){
      .argument_class(.data = Kinv, .class = "list")
      lapply(Kinv, function(Kinv) .check.objects.list(.data = Kinv, .class = c("matrix", "array", "ddiMatrix")))
    }
    if (!is.null(Q)){
      .argument_class(.data = Q, .class = c("matrix", "array", "ddiMatrix"))
    }
  }
  if (!user.mod){
    if (!is.null(Kinv)){
      lapply(Kinv, function(kinv){
        if (nrow(kinv) == ncol(kinv)){
          stop("At least one of the provided `Kinv` is not in the sparse form. Use `ASRgenomics::full2sparse()`.")
        }
        if (isFALSE(attributes(kinv)$INVERSE)){
          stop("The `attributes(Kinv)$INVERSE` must be `TRUE`.")
        }
      })
      condition <- lapply(Kinv, FUN = ASRgenomics:::Kinv.condition)
      ill <- which(unlist(condition) == "ill-conditioned")
      if (length(ill) > 0){
        stop(paste0("Kinv matrix(ces)", paste0(names(ill), collapse = " and "), " is/are 'ill-conditioned'."))
      }
      rm(condition, ill)
    }
  }
  if (message){
    message("Assuming `", names(pheno.data)[1], "` (first column of `pheno.data`) as genotype identifier.")
  }
  indiv <- names(pheno.data)[1]
  if (!user.mod){
    gen.no.vm <- gen[!gen %in% names(Kinv)]
    if (length(gen.no.vm) == 0 ) {gen.no.vm <- NULL}
    gen.vm <- gen[gen %in% names(Kinv)]
    if (length(gen.vm) == 0 ) {gen.vm <- NULL}
  }
  if (!user.mod){
    .variable_class(.data = pheno.data, .mandatory = TRUE, .variable = resp, .class = "numeric",
                          .class.action = "message", .mutate = TRUE, .message = message)
    mod.numeric <- c(weights, cov, total)
    .variable_class(.data = pheno.data, .mandatory = FALSE, .variable = mod.numeric, .class = "numeric",
                          .class.action = "message", .mutate = TRUE, .message = message)
    if (!is.null(fixedf)){
      fixedf_ <- unique(unlist(strsplit(x = fixedf, split = ":", fixed = TRUE)))
    } else {fixedf_ <- NULL}
    if (!is.null(randomf)){
      randomf_ <- unique(unlist(strsplit(x = randomf, split = ":", fixed = TRUE)))
    } else {randomf_ <- NULL}
    mod.factor <- c(gen, fixedf_, randomf_, residual)
    .variable_class(.data = pheno.data, .mandatory = FALSE, .variable = mod.factor,
                          .class = "factor", .class.action = "message",
                          .mutate = TRUE, .message = message)
    if(!is.null(residual)){
      if (message){
        message(paste0("Reordering `pheno.data` based on `", residual,
                       "` to use it as grouping factor for heterogeneous residual structure."))
      }
      ordered.pos <- order(pheno.data[[residual]])
      pheno.data <- pheno.data[ordered.pos, ]
    }
    rm(mod.numeric, mod.factor, fixedf_, randomf_)
  }
  if (!user.mod){
    if (any(!names(Kinv) %in% gen))
      stop("Names of objects in `Kinv` should correspond to values in `gen`.")
    if (any(duplicated(gen))){
      stop("Please, make sure that the objects in `gen` list are pointing to different variables in `pheno.data`.")
    }
    if (any(!is.na(match(gen, randomf)))){
      stop("Please, make sure to not pass the same unique factor to `gen` and `randomf`.")
    }
    if (length(grep(pattern = ":", x = gen, fixed = TRUE))){
      stop("Interactions are not allowed in `gen`; please, add viable interaction to `randomf`.")
    }
  }
  if (all(!is.null(weights), !is.null(residual))){
    stop("Use either `weights` or `residual` but not both.")
  }
  if (family == "binomial" & (!is.null(weights) | !is.null(residual))){
    stop("Residual `weights` or `residual` cannot be used for binomialy distributed data.")
  }
  if (!user.mod){
    if (caller == "gwas")
      na.prone.vars <- c(resp, weights)
    else
    if (caller == "gblup")
      na.prone.vars <- c(weights)
    if(is.data.table(pheno.data)){
      not.na.cases <- complete.cases(pheno.data[, ..na.prone.vars])
    } else {
      not.na.cases <- complete.cases(pheno.data[, na.prone.vars])
    }
    if (any(!not.na.cases)){
      pheno.data <- pheno.data[not.na.cases,]
      if (message){
        message("A total of ", sum(!not.na.cases),
          " sample(s) removed due to missing values in ",
          paste0(na.prone.vars, collapse = " or "), ".")
      }
    }
    if (caller == "gwas") {
      levels.gen.pre <- levels(pheno.data[[indiv]])
      pheno.data <- droplevels(pheno.data)
      levels.gen.pos <- levels(pheno.data[[indiv]])
      if (!identical(levels.gen.pre, levels.gen.pos)){
        stop("Genotypes dropped after 'NA' removal. Please deal with non-informative ",
          "levels before calling function 'gwas.asreml/gblup.asreml'.")
      }
    }
    rm(not.na.cases)
  }
  if (family == "binomial"){
    message("Checking binomial data:")
    if (!all(pheno.data[[resp]] %% 1 == 0)){
      stop("Only integers accepted in response varible if `family = 'binomial'`.")
    }
    if (any(pheno.data[[resp]] > 1)){
      if (is.null(total)){
        stop("The `total` counts of each sample must be provided when response varible is not binary and family is 'binomial'.")
      } else {
        message("  Response varible interpreted as number of successes in '", total, "' trials.")
      }
    } else {
      message("  Response varible interpreted as binary (0,1).")
      if (!is.null(total)){
        stop("Variable `total` must be `NULL` for binary data.")
      }
    }
  }
  match.data(
    pheno.data = pheno.data,
    indiv = indiv,
    map.data = map.data,
    geno.data = if(is.null(geno.data)) NULL else list(geno.data),
    Kinv = Kinv,
    caller = caller,
    Q = Q,
    message = message
  )
  if (!is.null(geno.data)){
    map.data$maf = ASRgenomics:::maf(geno.data)
  }
  if (bonferroni & !is.null(geno.data)){
    if (message) {
      message(col_blue(paste0('\nApplying Bonferroni correction to `pvalue.thr`.')))
    }
    pvalue.thr = 1 - (1 - pvalue.thr)^(1/ncol(geno.data))
  }
  if (!identical(as.character(pheno.data[[indiv]]), rownames(geno.data)) & !is.null(geno.data)){
    geno.data <- geno.data[as.character(pheno.data[[indiv]]), ]
  }
  if (!user.mod){
    if (!is.null(Q)) {
      Q <- Q[, 1:npc]
      colnames(Q) <- paste('Q', 1:npc, sep = '')
      pheno.data <- cbind.data.frame(pheno.data, Q[as.character(pheno.data[[indiv]]), ])
    } else {
      npc <- 0
    }
    if (npc > 0){
      group.cols <- c(ncol(pheno.data) - npc + 1, ncol(pheno.data))
    } else {
      group.cols <- NULL
    }
  }
  if (!user.mod){
    if (family == "gaussian"){
      pheno.data[[resp]] <- c(scale(pheno.data[[resp]], center = TRUE, scale = FALSE))
    }
  }
  if (!user.mod){
    code.asr <- call.generator(
      resp = resp, fixedf = fixedf, group.cols = group.cols, cov = cov,
      randomf = c(gen.no.vm, randomf), randomf.vm = gen.vm,
      vm.matrix.name = Kinv.name, data.name = "pheno.data",
      residual = residual, weights = weights,
      family = family, dispersion = dispersion, total = total,
      workspace = workspace)
  }
  if (!user.mod){
    if (message) {
      if (caller == "gwas"){
        message(col_blue('\nFitting base model.'))
      } else {
        message(col_blue('\nFitting model.'))
      }
    }
    if (P3D & family == "gaussian"){
      asreml::asreml.options(design = TRUE, trace = TRUE)
    } else {
      asreml::asreml.options(design = FALSE, trace = TRUE)
    }
    cur_log <- .try_asreml(
      x = "mod",
      call = code.asr,
      updates = 5,
      message = message
    )
    if (!is.null(cur_log$error)) {
      return(cur_log["error"])
    }
  }
  if ((!is.null(mod$call$group)) & (!user.mod)) {
    mod$call <- fix.call.grp_(mod_ = mod)
  }
  if (message){
    message(col_blue(paste0('Requesting Wald test.')))
  }
  .silence(aov <- tryCatch(
    expr = {asreml::wald(
      mod, denDF = "numeric", ssType = 'incremental', trace = FALSE)},
    error = function(msg) {return(NULL)}))
  if (is.null(aov) & message){
    message ("Unable to run Wald tests!")
  }
  if (!user.mod & (length(gen) > 0)){
    if (message){
      message(col_blue("\nObtaining heritability estimates."))
    }
    if(family == "gaussian"){
      skip.vc = !is.null(weights)
    }
    if (family == "binomial"){
      skip.vc = TRUE
    }
    heritability = repeatability.wrap(
      mod = mod,
      herit.numerators = c(
        get.vm.termvar(gen.vm, source = Kinv.name),
        gen.no.vm),
      skip.vc = skip.vc)
  }
  if (!is.null(geno.data)) {
    gwas.object <-
      gwas.asreml.core(
        mod = mod,
        weights = weights,
        map.data = map.data,
        geno.data = geno.data,
        Kinv = Kinv,
        pvalue.thr = pvalue.thr,
        P3D = P3D,
        user.mod = user.mod,
        maxiter.update = maxiter.update,
        inverse.update = inverse.update,
        threads = threads,
        obj.parallel = obj.parallel,
        message = message
      )
  }
  if (is.null(geno.data)) {
    gwas.object <- list(mod = mod)
  }
  gwas.object$call <- gwas.object$mod$call <- fix.call.param_(mod_ = mod)
  if (exists("heritability")) gwas.object <- append(gwas.object, heritability)
  if (exists("aov")) gwas.object <- append(gwas.object, list("aov" = aov))
  class(gwas.object) <- c("gwas.asreml", "list")
  return(gwas.object)
}
