#' @keywords internal
"_PACKAGE"
#' @aliases ASRgwas
#'
## usethis namespace: start
#' @import ASRgenomics
#' @import parallel
#' @import cli
#' @import ggplot2
#' @import ggrepel
#' @import data.table
#' @import stringr
#' @import ASRtools
#' @importFrom stats pnorm pt cor aggregate as.formula complete.cases median model.matrix na.exclude na.omit reformulate var
#' @importFrom methods formalArgs getFunction is
#' @importFrom utils combn packageVersion
#' @importFrom Rcpp sourceCpp
#' @useDynLib ASRgwas, .registration = TRUE
## usethis namespace: end
NULL
utils::globalVariables(
  c(
    "<<-",
    "..var_",
    "marker",
    "chrom",
    "pos",
    "chrom_pos",
    "..na.prone.vars",
    "nudge.x",
    "nudge.y",
    "x",
    "p.value",
    "Genotype",
    "Phenotype",
    "MAF",
    "obs"
  )
)
