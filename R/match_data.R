#' Check if objects to be used for GWAS model match
#'
#' @param pheno.data A data frame with all relevant columns (factors and covariates)
#' and phenotypic responses to be used (default = \code{NULL}).
#' @param indiv A character with name of the column in \code{pheno.data}
#' with the genotypes (default = \code{NULL}).
#' @param map.data A data frame with each marker's name, chromosome, and position.
#' Variable names \strong{must} be "marker", "chrom", and "pos" (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p},
#' with \eqn{n} individuals and \eqn{p} markers. Individual and marker names are
#' assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as: 0, 1, 2 (default = \code{NULL}).
#' @param Kinv A list of matrices representing the \strong{inverse} of the genomic relationship
#' matrix \eqn{\boldsymbol{Kinv}}. The matrices must be in sparse form (default = \code{NULL}).
#' @param Q A matrix of size \eqn{n} x \eqn{nv} with \eqn{nv} population structure-related covariates (\emph{e.g.},
#' principal component scores; default = \code{NULL}).
#' @param caller A character with the name of the function calling this function (default = \code{"gwas"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details
#' The matching performed on non-null objects are:
#' (i) between \'rownames' in \'geno.data' if \'geno.data' is a list;
#' (ii) between \'colnames' in \'geno.data' if \'geno.data' is a list;
#' (iii) between the \'rowNames' attribute of \'Kinv' if \'Kinv' is a list;
#' (iv) between the \'colNames' attribute of \'Kinv' if \'Kinv' is a list;
#' (v) between the levels of \'pheno.data[['indiv']]' and the \'rownames' of \'geno.data',
#' \'rowNames' and \'colNames' of \'Kinv, and the \'rownames' of \'Q'.
#' (vi) between the \'colnames' of \'geno.data' and \'map.data[['marker']]'.
#'
#' @return A hard stop if any pair of the provided objects do not match.
#'
#' @keywords internal
match.data <- function(
    pheno.data = NULL,
    indiv = NULL,
    map.data = NULL,
    geno.data = NULL,
    Kinv = NULL,
    Q = NULL,
    caller = c("gwas", "gblup"),
    message = TRUE) {
  if (!is.null(geno.data)) {
    if (length(geno.data) > 1) {
      match.in.list(container = geno.data, code = "rownames(obj)")
      match.in.list(container = geno.data, code = "colnames(obj)")
    }
  }
  if (!is.null(Kinv)) {
    if (length(Kinv) > 1) {
      match.in.list(container = Kinv, code = "attributes(obj)$rowNames")
      match.in.list(container = Kinv, code = "attributes(obj)$colNames")
    }
  }
  if (message)
    message("Using `levels(pheno.data[, indiv])` as reference for genotypes.")
  reference <- levels(pheno.data[[indiv]])
  if (!is.null(geno.data)) {
    in_pheno_not_geno <- reference %in% rownames(geno.data[[1]])
    in_geno_not_pheno <- rownames(geno.data[[1]]) %in% reference
    if (!all(in_pheno_not_geno))
      stop("`", indiv, "` has levels in `pheno.data` that are missing in `geno.data`")
    if (!all(in_geno_not_pheno))
      stop("`geno.data` has levels that are missing in `", indiv, "` of `pheno.data`")
  }
  if (!is.null(map.data) & !is.null(geno.data)) {
    if (!identical(map.data$marker, colnames(geno.data[[1]]))) {
      stop("`map.data$maker` and `colnames(geno.data)` are not identical.")
    }
  }
  if (!is.null(Kinv)) {
    in_pheno_not_Kinv <- reference %in% attributes(Kinv[[1]])$rowNames
    in_Kinv_not_pheno <- attributes(Kinv[[1]])$rowNames %in% reference
    if (!all(in_pheno_not_Kinv))
      stop("`", indiv, "` has levels in `pheno.data` that are missing in `Kinv`")
    if (!all(in_Kinv_not_pheno) && caller == "gblup" && message)
      message("`Kinv` has levels that are missing in `", indiv, "` of `pheno.data`")
    if (!all(in_Kinv_not_pheno) && caller == "gwas")
      stop("`Kinv` has levels that are missing in `", indiv, "` of `pheno.data`")
  }
  if (!is.null(Q)) {
    if (!identical(rownames(Q), reference)) {
      stop("`", indiv, "` has levels in `pheno.data` unmatched in `rownames(Q)`.")
    }
  }
  if (message) {
    message("No mismatches found in the provided data frames.")
  }
}
#' Match objects within a list (trap)
#'
#' @param list A list of objects to compare (default = \code{NULL}).
#' @param code The string of code to be used to check across objects (default = \code{NULL}).
#'
#' @details
#' The word "obj" is used in \code{code} to identify an
#' arbitrary object (each object in \code{list}).
#'
#' @return A hard stop if a mismatch is found between any pair of objects.
#'
#' @keywords internal
match.in.list <- function(container = NULL, code = NULL) {
  combinations <- combn(seq_along(container), m = 2, simplify = FALSE)
  matches <- lapply(combinations, function(X) {
    identical(
      eval(parse(
        text = gsub(x = code, pattern = "obj", replacement = "container[[X[1]]]")
      )),
      eval(parse(
        text = gsub(x = code, pattern = "obj", replacement = "container[[X[2]]]")
      ))
    )
  })
  if (!all(unlist(matches))) {
    stop("For at least one pair of \'obj' in ", deparse(substitute(list)), " list, the call \'", code, "' is not identical.")
  }
}
