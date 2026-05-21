#' Prepares and audits phenotypic and genomic data for GBLUP or GWAS
#'
#' @description This function checks the phenotypic and genomic data and to
#'   prepare the required data frames for downstream GBLUP or GWAS with
#'   \code{gblup.asreml} and \code{gwas.asreml}, respectively.
#'   The mandatory input is \code{pheno.data},
#'   \code{indiv}, and \code{geno.data}. If not provided, the map, and the
#'   matrices \eqn{\boldsymbol{K}}, \eqn{\boldsymbol{Kinv}} and
#'   \eqn{\boldsymbol{Q}} will be generated.
#'
#' @param pheno.data A mandatory data frame with all relevant columns (factors and
#'   covariates) and one or more phenotypic responses to be used (default =
#'   \code{NULL}).
#' @param indiv A mandatory character with the name of the column in \code{pheno.data}
#'   with the identification of treatments (genotypes) (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p}, with
#'   \eqn{n} individuals and \eqn{p} markers. Individual and marker names are
#'   assigned to \code{rownames} and \code{colnames}, respectively. SNP data is
#'   coded as: 0, 1, 2 (default = \code{NULL}).
#' @param resp A character with the name of the numeric variable in
#'   \code{pheno.data} with the response variable (default = \code{NULL}).
#' @param weights A character with the name of the numeric variable in
#'   \code{pheno.data} with the weights of each prediction (\emph{i.e.}, mean estimate)
#'   (default = \code{NULL}).
#' @param map.data An optional data frame with each marker's name, chromosome, and
#'   position. Variable names \strong{must} be: "marker", "chrom", and "pos"
#'   (default = \code{NULL}).
#' @param rename.markers If \code{TRUE} mathematical operators found on the
#'   markers' names will be replaced with the character "_" (default =
#'   \code{TRUE}).
#' @param K An optional symmetric matrix representing the genomic relationship
#'   or kinship matrix \eqn{\boldsymbol{K}}. The matrix must be in full form
#'   with size \eqn{n \times n}, with \eqn{n} individuals (genotypes). Genotypes
#'   names must be present as row and column names. If not
#'   provided, they will be obtained from \code{geno.data} (default =
#'   \code{NULL}).
#' @param Q.method A character with the data source to be used for the
#'   calculation of the \eqn{\boldsymbol{Q}} matrix by principal component
#'   analyses (PCA). Options are \code{"K"} (for kinship) and \code{"geno.data"}
#'   (for marker data) (default = \code{"K"}).
#' @param return A character vector specifying what to return. Options are:
#'   \code{"pheno.data"}, \code{"geno.data"}, \code{"K"}, \code{"Kinv"},
#'   \code{"map.data"}, \code{"Q"}, \code{"plot.scree"}, and
#'   \code{"eigenvalues"} (default = \code{"Kinv"}).
#' @param message If \code{TRUE}, diagnostic messages are printed on screen
#'   (default = \code{TRUE}).
#' @param ellipsis Further arguments to be called. Run \code{pre.gblup()} or
#'   \code{pre.gwas()} to check possible arguments and their default values.
#'   These arguments are called from functions:
#'   \link[ASRgenomics]{qc.filtering}, \link[ASRgenomics]{snp.pruning},
#'   \link[ASRgenomics]{G.matrix}, \link[ASRgenomics]{G.tuneup},
#'   \link[ASRgenomics]{kinship.diagnostics}, \link[ASRgenomics]{G.inverse},
#'   \link[ASRgenomics]{kinship.pca}, and \link[ASRgenomics]{snp.pca} available
#'   from the package \pkg{ASRgenomics}.
#'
#' @details
#'
#' A summary of the workflow considered in \code{pre.gblup} is:
#' \itemize{
#'   \item 1: Rename markers (if required).
#'   \item 2: Organize variables/columns in \code{pheno.data} (if necessary).
#'   \item 3: Sort provided data alphanumerically based on genotypes/individuals names.
#'   \item 4: Pre-match \code{pheno.data} and \code{geno.data} if the latter was provided.
#'   \item 5: Run \link[ASRgenomics]{qc.filtering} on \code{geno.data} to obtain filtered marker matrix.
#'   \item 6: Round genotypic values on \code{geno.data} (if decimals are present) as required by \link[ASRgenomics]{G.matrix}.
#'   \item 7: Obtain matrix \eqn{\boldsymbol{K}} using function \link[ASRgenomics]{G.matrix}, or use the one provided by user.
#'   \item 8: Run \link[ASRgenomics]{G.tuneup} on \eqn{\boldsymbol{K}} with blending, bending, and aligning (if requested by user).
#'   \item 9: Perform diagnostics on \eqn{\boldsymbol{K}} matrix using \link[ASRgenomics]{kinship.diagnostics}.
#'   \item 10: Repeat match of all datasets after \eqn{\boldsymbol{K}} quality control.
#'   \item 11: Reorder all matrices based on \code{levels(pheno.data[[indiv]])}.
#'   \item 12: Obtain inverse of kinship matrix (\eqn{\boldsymbol{Kinv}}) using \link[ASRgenomics]{G.inverse}.
#'   \item 13: If \eqn{\boldsymbol{Kinv}} is ill-conditioned (changes are cumulative) then:
#'     \itemize{
#'       \item 13.1: Run \link[ASRgenomics]{G.tuneup} on \eqn{\boldsymbol{K}} with blending of 0.05 using an identity matrix. If inverse is still ill-conditioned then go to 13.2; otherwise go to 14.
#'       \item 13.2: Run \link[ASRgenomics]{G.tuneup} on \eqn{\boldsymbol{K}} again with blending of 0.10 using an identity matrix. If inverse is still ill-conditioned then stop; otherwise go to 14.
#'     }
#'   \item 14: Return requested objects.
#' }
#'
#' A summary of the workflow considered in \code{pre.gwas} is:
#' \itemize{
#'   \item 1: Rename markers (if required).
#'   \item 2: Remove individual records with \code{NA} on \code{resp} and \code{weights}.
#'   \item 3: Organize variables/columns in \code{pheno.data} (if necessary).
#'   \item 4: Pre-match \code{pheno.data} and \code{geno.data}.
#'   \item 5: Run \link[ASRgenomics]{qc.filtering} on \code{geno.data} to obtain filtered marker matrix.
#'   \item 6: Round genotypic values on \code{geno.data} (if decimals are present) as required by \link[ASRgenomics]{G.matrix}.
#'   \item 7: Obtain matrix \eqn{\boldsymbol{K}} using function \link[ASRgenomics]{G.matrix}, or use the one provided by user.
#'   \item 8: Run \link[ASRgenomics]{G.tuneup} on \eqn{\boldsymbol{K}} with blending, bending, and aligning (if requested by user).
#'   \item 9: Perform diagnostics on \eqn{\boldsymbol{K}} matrix using \link[ASRgenomics]{kinship.diagnostics}.
#'   \item 10: Repeat match of all datasets after \eqn{\boldsymbol{K}} quality control.
#'   \item 11: Reorder all matrices based on \code{levels(pheno.data[[indiv]])}.
#'   \item 12: Obtain inverse of kinship matrix (\eqn{\boldsymbol{Kinv}}) using \link[ASRgenomics]{G.inverse}.
#'   \item 13: If \eqn{\boldsymbol{Kinv}} is ill-conditioned (changes are cumulative) then:
#'     \itemize{
#'       \item 13.1: Run \link[ASRgenomics]{G.tuneup} on \eqn{\boldsymbol{K}} with blending of 0.05 using an identity matrix. If inverse is still ill-conditioned then go to 13.2; otherwise go to 14.
#'       \item 13.2: Run \link[ASRgenomics]{G.tuneup} on \eqn{\boldsymbol{K}} again with blending of 0.10 using an identity matrix. If inverse is still ill-conditioned then stop; otherwise go to 14.
#'     }
#'   \item 14: Obtain matrix \eqn{\boldsymbol{Q}} based on \code{geno.data} or final \eqn{\boldsymbol{K}} matrix.
#' }
#'
#' The majority of the default values used for arguments are those suggested
#' in \pkg{ASRgenomics}, but these can be modified via the ellipsis argument (\code{...}).
#' Nevertheless, some defaults have been changed for better implementation of
#' GBLUP and GWAS downstream procedures.
#'
#' The modified defaults for GBLUP are: \code{ncp = 20}, \code{maf = 0.05},
#' \code{marker.callrate = 0.2}, \code{ind.callrate = 0.2}, \code{clean.diagonal
#' = TRUE}, \code{clean.duplicate = TRUE}, \code{diagonal.thr.large = 1.8}, and
#' \code{diagonal.thr.small = 0.2}. It might be important to adjust these
#' parameters along with GBLUP runs to identify the best analytical approach.
#'
#' The modified defaults for GWAS are: \code{ncp = 20}, \code{maf = 0.05},
#' \code{marker.callrate = 0.2}, \code{ind.callrate = 0.2}, \code{clean.diagonal
#' = TRUE}, \code{clean.duplicate = TRUE}, \code{diagonal.thr.large = 1.8}, and
#' \code{diagonal.thr.small = 0.2}. It might be important to adjust these
#' parameters along with GWAS runs to identify best analytical approach.
#'
#' For GWAS, the \eqn{\boldsymbol{Q}} matrix can be calculated based on
#' \eqn{\boldsymbol{K}} or \code{geno.data}. Nevertheless, if missing values are
#' still present in \code{geno.data} after cleaning, \code{Q.method = "K"} will
#' be always used, as this method allows for missing values.
#'
#' For more information on the aforementioned functions, please refer to the
#' \pkg{ASRgenomics} help.
#'
#' @return A list with several data frames and objects to be used in downstream
#'   GBLUP and GWAS model fits.
#'
#' For GBLUP, these items are:
#'
#' \itemize{
#'   \item \code{pheno.data}: A verified phenotypic data frame with all relevant columns (and rows) and responses to be used.
#'   \item \code{map.data}: A verified (or generated) data frame with markers' names, chromosomes, and positions.
#'   \item \code{geno.data}: A filtered and verified matrix with marker data.
#'   \item \code{K}: A list of verified and tuned-up kinship matrices.
#'   \item \code{Kinv}: A list of verified and tuned-up inverse kinship matrices in sparse format (ready for ASReml-R).
#' }
#'
#' For GWAS, these items are:
#'
#' \itemize{
#'   \item \code{pheno.data}: A verified phenotypic data frame with all relevant columns (and rows) and responses to be used.
#'   \item \code{map.data}: A verified (or generated) data frame with markers' names, chromosomes, and positions.
#'   \item \code{geno.data}: A filtered and verified matrix with marker data.
#'   \item \code{K}: A verified and tuned-up kinship matrix (if requested).
#'   \item \code{Kinv}: A verified and tuned-up inverse of the kinship matrix in sparse format (ready for ASReml-R).
#'   \item \code{Q}: A population structure \eqn{\boldsymbol{Q}} matrix of principal component scores.
#'   \item \code{plot.scree}: The scree plot for the components in the \eqn{\boldsymbol{Q}} matrix.
#'   \item \code{eigenvalues}: The eigenvalues for the components in the \eqn{\boldsymbol{Q}} matrix.
#' }
#'
#' @keywords internal
pre.generic <- function(
    pheno.data = NULL,
    indiv = NULL,
    resp = NULL,
    weights = NULL,
    geno.data = NULL,
    map.data = NULL,
    rename.markers = TRUE,
    K = NULL,
    Q.method = c("none", "K", "geno.data"),
    return = c("pheno.data", "Kinv"),
    message = TRUE,
    caller = c("gwas", "gblup"),
    ellipsis = NULL
) {
  caller <- match.arg(caller)
  if (is.null(indiv) & is.null(geno.data) &
      is.null(map.data) &
      is.null(pheno.data)) {
    return(asrgen.args())
  }
  Q.method <- match.arg(Q.method)
  return.K <- ifelse({"K" %in% return}, yes = TRUE, no =  FALSE)
  return.Kinv <- ifelse({"Kinv" %in% return}, yes = TRUE, no =  FALSE)
  return.geno.data <- ifelse({"geno.data" %in% return}, yes = TRUE, no =  FALSE)
  return.map.data <- ifelse({"map.data" %in% return}, yes = TRUE, no =  FALSE)
  return.Q <- ifelse({"Q" %in% return}, yes = TRUE, no =  FALSE)
  calculate.K <- ifelse(return.K || return.Kinv || (Q.method ==   "K"),
                        yes = TRUE, no = FALSE)
  if (is.null(geno.data) & return.geno.data)
    stop("Argument `geno.data` is mandatory when its return is requested.",
         call. = FALSE)
  if (message) {message(col_blue("\nChecking input."))}
  .argument_class(.data = pheno.data, .class = "data.frame")
  pheno.data <- droplevels(pheno.data)
  .variable_class(
    .data = pheno.data, .mandatory =  TRUE, .variable = indiv, .class = "factor",
    .class.action = "message", .mutate= TRUE, .message = message)
  .variable_class(
    .data = pheno.data, .mandatory =  FALSE, .variable = resp, .class = "numeric",
    .class.action = "message", .mutate= TRUE, .message = message)
  .variable_class(
    .data = pheno.data, .mandatory =  FALSE, .variable = weights, .class = "numeric",
    .class.action = "message", .mutate= TRUE, .message = message)
  if (!is.null(geno.data)) {
    .argument_class(.data = geno.data, .class = c("matrix", "array"))
    if(is.null(rownames(geno.data))){
      stop("Individual names not assigned to rows of `geno.data`.")
    }
    if(is.null(colnames(geno.data))){
      stop("Marker names not assigned to columns of `geno.data`.")
    }
  }
  .argument_class(.data = message, .class = "logical")
  if (!is.null(map.data)) {
    .argument_class(.data = map.data, .class = "data.frame")
    if (!identical(names(map.data), c("marker", "chrom", "pos")))
      stop("Column names on `map.data` must be: 'marker', 'chrom', and 'pos'.")
  }
  if (!is.null(K)) {
    .argument_class(.data = K, .class = c("matrix", "array"))
    if (is.null(rownames(K)) |
        is.null(colnames(K)) |
        !identical(rownames(K), colnames(K))) {
      stop("Individual names not assigned to rows and columns of matrix K.")
    }
  }
  if (rename.markers & !is.null(geno.data)) {
    fixed_marker_names <- .replace_operators(colnames(geno.data))
    if (!is.null(fixed_marker_names))
      colnames(geno.data) <- fixed_marker_names
    if (!is.null(map.data)) {
      fixed_marker_names <- .replace_operators(map.data$marker)
      if (!is.null(fixed_marker_names))
        map.data$marker <- fixed_marker_names
    }
    if (exists("fixed_marker_names", inherits = FALSE)) rm(fixed_marker_names)
  }
  if (message) {
    message(col_blue("\nChecking data frame `pheno.data`."))
  }
  if (caller == "gwas")
    na.prone.vars <- c(resp, weights)
  if (caller == "gblup")
    na.prone.vars <- NULL
  if (!is.null(na.prone.vars)){
    if(is.data.table(pheno.data)){
      not.na.cases <- complete.cases(pheno.data[, ..na.prone.vars])
    } else {
      not.na.cases <- complete.cases(pheno.data[, na.prone.vars])
    }
    if (any(!not.na.cases)){
      pheno.data <- pheno.data[not.na.cases,]
      if (message){
        message("A total of ", sum(!not.na.cases),
                " sample(s) were removed due to missing values in ",
                paste0(na.prone.vars, collapse = " or "), ".")
      }
    }
    rm(not.na.cases)
  }
  pheno.data <- droplevels(pheno.data)
  if(match(indiv, names(pheno.data)) != 1){
    if (message) {
      message("Reordering columns in `pheno.data` so `indiv` ",
              "is the first column [required by `gwas.asreml()`].")
    }
    is.indiv <- names(pheno.data) %in% indiv
    pheno.data <- pheno.data[, c(which(is.indiv), which(!is.indiv))]
  }
  if (message) {
    message(col_blue("\nSorting datasets if needed."))
  }
  ordered_levels <- str_sort(pheno.data[[indiv]], numeric = TRUE) |> na.omit() |>
    unique()
  if (!identical(levels(pheno.data[[indiv]]), ordered_levels)) {
    if (message)
      message("Reordering levels in `pheno.data[, indiv]` alphanumerically.")
    pheno.data[[indiv]] <- factor(pheno.data[[indiv]], levels = ordered_levels)
  }
  if (!is.null(geno.data)) {
    ordered_levels <- str_sort(rownames(geno.data), numeric = TRUE) |>
      na.omit()
    if (!identical(rownames(geno.data), ordered_levels)) {
      message("Reordering rows in `geno.data` alphanumerically.")
      geno.data <- geno.data[ordered_levels, , drop = FALSE]
    }
  }
  if (exists("ordered_levels", inherits = FALSE))
    rm(ordered_levels)
  if (!is.null(K)) {
    ordered_levels <- str_sort(rownames(K), numeric = TRUE)
    if (!identical(rownames(K), ordered_levels)) {
      message("Reordering rows and columns in `K` alphanumerically.")
      K <- K[ordered_levels, ordered_levels]
    }
  }
  if (!is.null(geno.data)) {
    if (message) {
      message(col_blue("\nMatching individuals in `pheno.data` ",
                       "and `geno.data`."))
    }
    indiv.names <- levels(pheno.data[[indiv]])
    unique.phen <- !indiv.names %in% rownames(geno.data)
    unique.geno.data <- !rownames(geno.data) %in% indiv.names
    if (any(c(unique.phen, unique.geno.data))){
      if (message){
        message("Removing ", sum(unique.phen),
                " genotypes present in `pheno.data` and absent in ",
                "`geno.data`.")
        if (caller == "gwas")
          message("Removing ", sum(unique.geno.data),
                  " genotypes present in `geno.data` and absent ",
                  "in `pheno.data`.")
        if (caller == "gblup")
          message("Found ", sum(unique.geno.data),
                  " genotypes present in `geno.data` and absent ",
                  "in `pheno.data`.")
      }
      pheno.data <- pheno.data[pheno.data[[indiv]] %in%
                                 indiv.names[!unique.phen], ]
      pheno.data <- droplevels(pheno.data)
      if (caller == "gwas")
        geno.data <- geno.data[levels(pheno.data[[indiv]]), ]
    }
    if (!is.null(map.data)) {
      ord_by_chrom_n_pos <- with(map.data, order(chrom, pos))
      if (!identical(ord_by_chrom_n_pos, seq_len(nrow(map.data)))) {
        if (message)
          message("Ordering `geno.data` and `map.data` by `chrom` and `pos`.")
        map.data <- map.data[ord_by_chrom_n_pos, ]
        geno.data <- geno.data[, ord_by_chrom_n_pos]
      }
    }
    if (!is.null(K)) {
      if (!all(rownames(K) %in% rownames(geno.data)) ||
          !all(colnames(K) %in% rownames(geno.data)))
        stop("Please provide the `geno.data` matrix used to calculate `K`.")
    }
    rm(indiv.names, unique.phen, unique.geno.data)
  }
  qc.filtering <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::qc.filtering)
  G.matrix <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::G.matrix)
  G.tuneup <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::G.tuneup)
  kinship.diagnostics <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::kinship.diagnostics)
  G.inverse <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::G.inverse)
  kinship.pca <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::kinship.pca)
  snp.pca <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::snp.pca)
  snp.pruning <-
    ellipsis.subset(ellipsis = ellipsis, fun = ASRgenomics::snp.pruning)
  if (!is.null(geno.data)){
    if (message) {
      message(col_blue("\nQuality control of data frame `geno.data` ",
                       "in progress."))
    }
    if (any(dim(geno.data) == 0))
      stop("Empty `geno.data` matrix, possibly a marker or genotype name mismatch with other datasets.", call. = FALSE)
    if (!is.null(map.data)) {
      if (any(dim(map.data) == 0)) {
      stop("Empty `map.data` matrix, possibly a marker or genotype name mismatch with other datasets.", call. = FALSE)
      }
    }
    tmp.qc <- ASRgenomics::qc.filtering(
      M = geno.data,
      map = map.data,
      marker = "marker",
      chrom = "chrom",
      pos = "pos",
      base = qc.filtering$base,
      ref = qc.filtering$ref,
      marker.callrate = qc.filtering$marker.callrate,
      ind.callrate = qc.filtering$ind.callrate,
      maf = qc.filtering$maf,
      heterozygosity = qc.filtering$heterozygosity,
      Fis = qc.filtering$Fis,
      na.string = qc.filtering$na.string,
      impute = qc.filtering$impute,
      message = message,
      Mrecode = FALSE,
      plots = FALSE,
      digits = 2
    )
    geno.data <- tmp.qc$M.clean
    if (!is.null(map.data)){
      map.data <- tmp.qc$map
    }
    rm(tmp.qc)
  }
  if (!is.null(geno.data)){
    if (!isFALSE(ellipsis$pruning.thr)){
      geno.data <- ASRgenomics::snp.pruning(
        M = geno.data,
        map = map.data,
        marker = "marker",
        chrom = "chrom",
        pos = "pos",
        method = "correlation",
        criteria = snp.pruning$criteria,
        pruning.thr = snp.pruning$pruning.thr,
        by.chrom = snp.pruning$by.chrom,
        window.n = snp.pruning$window.n,
        overlap.n = snp.pruning$overlap.n,
        iterations = snp.pruning$iterations,
        seed = snp.pruning$seed,
        message = message
      )$Mpruned
    }
  }
  if (!is.null(geno.data)){
    if (message) {
      message(col_blue("\nChecking data frame `geno.data` for missing data."))
    }
    any.miss.geno.data <- any(is.na(geno.data))
    if (message & any.miss.geno.data) {
      message("Missing data present in `geno.data`. This will impact in ",
              "processing time of downstream analyses.")
      message("Use `impute = TRUE` to request mean-based imputation.")
    }
  }
  if (!is.null(geno.data)){
    if (message) {
      message(col_blue("\nMatching individuals in `pheno.data` and ",
                       "`geno.data`."))
    }
    indiv.names <- levels(pheno.data[[indiv]])
    unique.phen <- !indiv.names %in% rownames(geno.data)
    unique.geno.data <- !rownames(geno.data) %in% indiv.names
    if (any(c(unique.phen, unique.geno.data))){
      if (message){
        message("Removing ", sum(unique.phen),
                " genotypes present in `pheno.data` and absent in ",
                "\'geno.data'.")
        if (caller == "gwas")
          message("Removing ", sum(unique.geno.data),
                  " genotypes present in `geno.data` and absent ",
                  "in `pheno.data`.")
        if (caller == "gblup")
          message("Found ", sum(unique.geno.data),
                  " genotypes present in `geno.data` and absent ",
                  "in `pheno.data`.")
      }
      pheno.data <- pheno.data[pheno.data[[indiv]] %in%
                                 indiv.names[!unique.phen], ]
      pheno.data <- droplevels(pheno.data)
      if (caller == "gwas")
        geno.data <- geno.data[levels(pheno.data[[indiv]]), ]
    }
    if (!is.null(map.data)) {
      map.data <- map.data[map.data[["marker"]] %in% colnames(geno.data), ]
    }
    if (!is.null(K)) {
      sorting_crit <- match(rownames(geno.data), colnames(K))
      K <- K[sorting_crit, sorting_crit]
    }
    rm(indiv.names, unique.phen, unique.geno.data)
  }
  if (!is.null(geno.data)){
    if (message) {
      message(col_blue("\nObtaining relationship matrix (K)."))
    }
    if (!all(na.omit(c(geno.data) %% 1 == 0)))
    {
      if (message){
        message("Converting values in `geno.data` to integers.")
      }
      geno.data <- round(x = geno.data, digits = 0)
    }
  }
  if (!is.null(K)) {
    message("Using provided `K` matrix.")
  }
  else if (calculate.K) {
    if (message) {
      message(paste0("Computing K matrix with ", G.matrix$method), "'s method.")
    }
    if (is.null(geno.data)) stop(
      "K, Kinv, or K-based Q return requested but no geno.data or K provided!",
      call. = FALSE)
    .silence(
      .code =
        K <- ASRgenomics::G.matrix(
          M = geno.data,
          method = G.matrix$method,
          na.string = NA,
          sparseform = FALSE,
          digits = 8
        )$G
    )
  } else {
    K <- NULL
  }
  if ((G.tuneup$blend | G.tuneup$bend | G.tuneup$align) & !is.null(K)) {
    if (message) message(col_blue("\nImplementing tune-up on K matrix."))
    if (!is.null(G.tuneup$A)){
      tmpGA <- ASRgenomics::match.G2A(
        A = G.tuneup$A,
        G = K,
        clean = TRUE,
        ord = TRUE,
        mism = FALSE)
      G.tuneup$A <- tmpGA$Aclean
      K <- tmpGA$Gclean
    }
    K <- ASRgenomics::G.tuneup(
      G = K,
      A = G.tuneup$A,
      blend = G.tuneup$blend,
      pblend = G.tuneup$pblend,
      bend = G.tuneup$bend,
      eig.tol = G.tuneup$eig.tol,
      align = G.tuneup$align,
      rcn = FALSE,
      digits = 8,
      sparseform = FALSE,
      determinant = FALSE,
      message = message
    )$Gb
  }
  if (!is.null(K)) {
    if (message) message(col_blue("\nPerforming diagnostics on K matrix."))
    tmp.K <- ASRgenomics::kinship.diagnostics(
      K = K,
      diagonal.thr.large = kinship.diagnostics$diagonal.thr.large,
      diagonal.thr.small = kinship.diagnostics$diagonal.thr.small,
      duplicate.thr = kinship.diagnostics$duplicate.thr,
      clean.diagonal = TRUE,
      clean.duplicate = TRUE,
      message = message
    )
    if (!is.null(tmp.K$list.diagonal) && !kinship.diagnostics$clean.diagonal) {
      cat("\n")
      stop("\nOut of range diagonal values on K matrix, consider:\n",
           " > Increasing the range of allowed values (not ideal); or\n",
           " > Including `clean.diagonal = TRUE` to this function.",
           call. = FALSE)
    }
    if (!is.null(tmp.K$list.duplicate) && !kinship.diagnostics$clean.duplicate) {
      cat("\n")
      stop("\nDuplicates identified but not removed consider:\n",
           " > Increasing the threshold; or\n",
           " > Including `clean.duplicates = TRUE` to this function.",
           call. = FALSE)
    }
    if (!is.null(tmp.K$clean.kinship)){
      if (length(tmp.K$clean.kinship) <= 1){
        cat("\n")
        stop("No genotypes remain after cleaning based on diagonal values ",
             "and/or duplicates.", call. = FALSE)
      }
      if (message) message("Cleaning K matrix based on diagonal values ",
                           "and/or duplicates.")
      K <- tmp.K$clean.kinship
    }
    if (max(diag(K), na.rm = TRUE) < 0.99)
      warning("All diagonal values of K are lower than 0.99. This may indicate ",
              "problems.", call. = FALSE)
    rm(tmp.K)
  }
  if (!is.null(K)){
    if (message) {
      message(col_blue("\nMatching individuals in `pheno.data` and ",
                       "`K`."))
    }
    indiv.names <- levels(pheno.data[[indiv]])
    unique.phen <- !indiv.names %in% rownames(K)
    unique.geno.data <- !rownames(K) %in% indiv.names
    if (any(c(unique.phen, unique.geno.data))) {
      if (message){
        message("Removing ", sum(unique.phen),
                " genotypes present in `pheno.data` and absent in ",
                "`K`.")
        if (caller == "gwas")
          message("Removing ", sum(unique.geno.data),
                  " genotypes present in `K` and absent ",
                  "in `pheno.data`.")
        if (caller == "gblup")
          message("Found ", sum(unique.geno.data),
                  " genotypes present in `K` and absent ",
                  "in `pheno.data`.")
      }
    }
    pheno.data <- pheno.data[pheno.data[[indiv]] %in%
                               indiv.names[!unique.phen], ]
    pheno.data <- droplevels(pheno.data)
    indiv.names <- levels(pheno.data[[indiv]])
    if (!is.null(geno.data)) {
      geno.data <- geno.data[rownames(K), ]
      if (caller == "gwas")
        geno.data <- geno.data[indiv.names, ]
      if (!is.null(map.data))
        map.data <- map.data[map.data[["marker"]] %in% colnames(geno.data), ]
    }
    if (caller == "gwas")
      K <- K[indiv.names, indiv.names]
    rm(indiv.names, unique.phen, unique.geno.data)
  }
  if (return.Kinv) {
    Kinv.status <- ""
    if (message) message(col_blue("\nObtaining inverse of K matrix."))
    tryCatch(expr =
               Kinv <- ASRgenomics::G.inverse(
                 G = K,
                 rcn.thr = G.inverse$rcn.thr,
                 eig.tol = G.inverse$eig.tol,
                 digits = 8,
                 sparseform = FALSE,
                 message = message
               )$Ginv,
             error = function(holder)
               return(Kinv.status <<- "ill-conditioned"))
    if (exists("Kinv", inherits = FALSE)) {Kinv.status <- ASRgenomics:::Kinv.condition(Kinv)}
    if (Kinv.status == "ill-conditioned"){
      Kinv.status <- ""
      if (message) message(
        col_yellow(
          "The inverse of the K matrix is ","ill-conditioned.",
          " Implementing blending with p = 0.05 using an identity matrix."))
      K <- ASRgenomics::G.tuneup(
        G = K,
        A = NULL,
        blend = TRUE,
        pblend = .05,
        bend = FALSE,
        eig.tol = G.tuneup$eig.tol,
        align = FALSE,
        rcn = FALSE,
        digits = 8,
        sparseform = FALSE,
        determinant = FALSE,
        message = FALSE
      )$Gb
      tryCatch(expr =
                 Kinv <- ASRgenomics::G.inverse(
                   G = K,
                   eig.tol = G.inverse$eig.tol,
                   rcn.thr = G.inverse$rcn.thr,
                   digits = 8,
                   sparseform = FALSE,
                   message = message
                 )$Ginv,
               error = function(holder)
                 return(Kinv.status <<- "ill-conditioned"))
    }
    if (exists("Kinv", inherits = FALSE)) {Kinv.status <- ASRgenomics:::Kinv.condition(Kinv)}
    if (Kinv.status == "ill-conditioned"){
      Kinv.status <- ""
      if (message) message(col_yellow(
        "The inverse of the K matrix is still ill-conditioned.",
        " Implementing blending with p = 0.10 using an identity matrix."))
      K <- ASRgenomics::G.tuneup(
        G = K,
        A = NULL,
        blend = TRUE,
        pblend = .1,
        bend = FALSE,
        eig.tol = G.tuneup$eig.tol,
        align = FALSE,
        rcn = FALSE,
        digits = 8,
        sparseform = FALSE,
        determinant = FALSE,
        message = FALSE
      )$Gb
      tryCatch(expr =
                 Kinv <- ASRgenomics::G.inverse(
                   G = K,
                   rcn.thr = G.inverse$rcn.thr,
                   eig.tol = G.inverse$eig.tol,
                   digits = 8,
                   sparseform = FALSE,
                   message = message
                 )$Ginv,
               error = function(holder)
                 return(Kinv.status <<- "ill-conditioned"))
    }
    if (Kinv.status == "ill-conditioned")
      stop("\n\'Inverse of K matrix is still ill-conditioned after blending ",
           "with 0.05 and 0.10. User action is required.", call. = FALSE)
    Kinv <- list(ASRgenomics::full2sparse(Kinv))
    names(Kinv) <- indiv
  } else {
    Kinv <- NULL
  }
  if (return.Q) {
    if(Q.method == "K" || any.miss.geno.data){
      if (message){
        message(col_blue("\nObtaining Q matrix from K matrix."))
      }
      if (Q.method == "geno.data" & any.miss.geno.data){
        message("Obtaining Q matrix from K matrix as `geno.data` has ",
                "missing values.")
      }
      Q <- ASRgenomics::kinship.pca(
        K = K,
        scale = kinship.pca$scale,
        label = kinship.pca$label,
        ncp = kinship.pca$ncp,
        groups = NULL
      )
    } else {
      if (message){
        message(col_blue("\nObtaining Q matrix from `geno.data` matrix."))
      }
      Q <- ASRgenomics::snp.pca(
        M = geno.data,
        label = snp.pca$label,
        ncp = snp.pca$ncp,
        groups = NULL
      )
    }
  } else {
    Q = NULL
  }
  if (is.null(map.data) & !is.null(geno.data) & return.map.data) {
    if (message) message(col_blue("\nCreating dummy map."))
    map.data <- ASRgenomics:::dummy.map_(
      marker.id = colnames(geno.data),
      message = FALSE)
  }
  if (!return.geno.data) geno.data <- NULL
  if (!return.map.data) map.data <- NULL
  if (return.K){
    K <- list(K)
    names(K) <- indiv
  } else {
    K <- NULL
  }
  if (message) message(col_blue("\nReturning objects."))
  structure(
    list(
      pheno.data = pheno.data,
      geno.data = geno.data,
      map.data = map.data,
      Kinv = Kinv,
      K = K,
      Q = Q$pca.scores,
      plot.scree = Q$plot.scree,
      eigenvalues = Q$eigenvalues),
    class = "predata")
}
