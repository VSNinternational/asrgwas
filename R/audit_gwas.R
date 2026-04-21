#' Checks and audits phenotypic and genomic datasets prior to fitting GWAS models
#'
#' @description Implements checks and basic exploratory data analysis on
#' phenotypic and genomic datasets in order to ensure they will work adequately on
#' downstream GWAS model fits using the function \code{gwas.asreml}. Mismatches
#' within and between datasets are reported to help detect inconsistencies.
#'
#' This function can be used in place of \code{pre.gwas} in order to provide
#' \code{gwas.asreml} with all the required elements. These can be constructed
#' in other routines, libraries or your own code. This function will ensure that
#' all attributes and elements are suitable for use with \code{gwas.asreml}.
#'
#' @param pheno.data A mandatory data frame with all relevant columns (phenotypic
#'   response, factors and covariates) to be used with \code{gwas.asreml} (default =
#'   \code{NULL}).
#' @param indiv A mandatory character string of the name of the column in
#'   \code{pheno.data} with the identification of treatments (genotypes)
#'   (default = \code{NULL}).
#' @param resp A character string of the name of the numeric variable in
#'   \code{pheno.data} that identifies the response variable (default = \code{NULL}).
#' @param Kinv A list of matrices representing the \strong{inverse} of
#'   the genomic relationship matrix \eqn{\boldsymbol{Kinv}}. The matrices must
#'   be in sparse form (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p},
#'   with \eqn{n} individuals and \eqn{p} markers. Individual and marker names
#'   are assigned to \code{rownames} and \code{colnames}, respectively. SNP data
#'   is coded as: 0, 1, 2, and \code{NA} (default = \code{NULL}).
#' @param map.data An optional data frame with the marker names, chromosomes,
#'   and positions. Variable names \strong{must} be: "marker", "chrom", and "pos"
#'   (default = \code{NULL}).
#' @param Q An optional matrix of size \eqn{n} x \eqn{nv}, with \eqn{nv}
#'   population structure-related covariates (\emph{e.g.}, principal component scores;
#'   default = \code{NULL}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen
#'   (default = \code{TRUE}).
#'
#' @return A list containing information about the dataset, and summary
#'   statistics for the phenotypic, relationship and molecular data.
#'
#' \itemize{
#'   \item \code{trial.stats}: Summary statistics for the response variable.
#'   \item \code{gen.stats}: Summary statistics for each of the genotypes.
#'   \item \code{n.genotypes}: Total number of genotypes to be used downstream.
#'   \item \code{n.markers}: Total number of markrs to be used downstream.
#'   \item \code{ind.callrate.stats}: Minimum and maximum values for call rates for individuals.
#'   \item \code{marker.callrate.stats}: Minimum and maximum values for call rates for markers.
#'   \item \code{maf.stats}: Minimum and maximum values for minor allele frequency (MAF).
#'   \item \code{Kinv.condition}: Informs on the condition of the \eqn{\boldsymbol{Kinv}} matrix.
#' }
#'
#' @details The number of genotypes is tentatively obtained from the provided
#' input with the following priority order: \code{geno.data}, \code{Kinv},
#' levels of \code{pheno.data[[indiv]]}.
#'
#' @export
#'
#' @examples
#'
#' # Prepare Apricot datasets (simple run).
#' gwas.data <- pre.gwas(
#'   pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
#'   geno.data = geno.apricot)
#'
#' # Audit data frames and matrices prior to GWAS.
#' audit <- audit.gwas(
#'  pheno.data = gwas.data$pheno.data, resp = "Sucrose", indiv = "Ind",
#'  Kinv = gwas.data$Kinv, geno.data = gwas.data$geno.data, Q = gwas.data$Q)
#' audit$trial.stats
#' audit$gen.stats
#' audit$n.genotypes
#' audit$n.markers
#' audit$ind.callrate.stats
#' audit$marker.callrate.stats
#' audit$maf.stats
#' audit$Kinv.condition
audit.gwas <- function(
    pheno.data = NULL,
    indiv = NULL,
    resp = NULL,
    Kinv = NULL,
    geno.data = NULL,
    map.data = NULL,
    Q = NULL,
    message = TRUE){
  audit.call <- sys.call()
  audit.call[[1]] <- str2lang("audit.generic")
  audit.call$caller <- "gwas"
  eval(audit.call)
}
