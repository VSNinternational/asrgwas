#' Checks and prepares phenotypic and genomic datasets prior to fitting GWAS models
#'
#' @description This function checks, corrects and matches the phenotypic and genomic data to
#' prepare the required data frames for downstream GWAS with \code{gwas.asreml}.
#' If not provided, the map, and the matrices \eqn{\boldsymbol{K}},
#' \eqn{\boldsymbol{Kinv}} and \eqn{\boldsymbol{Q}} will be automatically generated.
#'
#' @param pheno.data A mandatory data frame with all relevant columns (phenotypic
#'   response, factors and covariates) to be used with \code{gwas.asreml} (default =
#'   \code{NULL}).
#' @param indiv A mandatory character string of the name of the column in \code{pheno.data}
#'   that identifies treatments (genotypes) (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p}, with
#'   \eqn{n} individuals and \eqn{p} markers. Individual and marker names are
#'   assigned to \code{rownames} and \code{colnames}, respectively. SNP data is
#'   coded as: 0, 1, 2, and \code{NA} (default = \code{NULL}).
#' @param resp A character string of the name of the numeric variable in
#'   \code{pheno.data} that identifies the response variable (default = \code{NULL}).
#' @param weights A character string of the name of the numeric variable in
#'   \code{pheno.data} that specifies the weights of each prediction (\emph{i.e.}, mean estimate)
#'   (default = \code{NULL}).
#' @param map.data A data frame with the marker names, chromosomes, and
#'   positions. Variable names \strong{must} be: "marker", "chrom", and "pos"
#'   (default = \code{NULL}).
#' @param K A symmetric matrix representing the genomic relationship
#'   or kinship matrix \eqn{\boldsymbol{K}}. The matrices must be in full form
#'   with size \eqn{n \times n} with \eqn{n} individuals (genotypes). Genotypes
#'   names must be present as row and column names. If not
#'   provided, they will be obtained from \code{geno.data} (default =
#'   \code{NULL}).
#' @param Q.method A character with the data source to be used for the
#'   calculation of the \eqn{\boldsymbol{Q}} matrix by principal component
#'   analyses (PCA). Options are \code{"K"} (for kinship) and \code{"geno.data"}
#'   (for marker data) (default = \code{"K"}).
#' @param rename.markers If \code{TRUE} mathematical operators found on the
#'   markers' names will be replaced with the character "_" (default =
#'   \code{TRUE}).
#' @param return A character vector specifying what to return. Options are:
#'   \code{"pheno.data"}, \code{"geno.data"}, \code{"K"}, \code{"Kinv"},
#'   \code{"map.data"}, \code{"Q"}, \code{"plot.scree"}, and
#'   \code{"eigenvalues"} (default = all).
#' @param message If \code{TRUE}, diagnostic messages are printed on screen
#'   (default = \code{TRUE}).
#' @param ... Further arguments to be called from \pkg{ASRgenomics}. Run the command
#'   \code{"pre.gwas()"} for a list of all possible arguments. These arguments are
#'   called from the following functions:
#'   \link[ASRgenomics]{qc.filtering}, \link[ASRgenomics]{snp.pruning},
#'   \link[ASRgenomics]{G.matrix}, \link[ASRgenomics]{G.tuneup},
#'   \link[ASRgenomics]{kinship.diagnostics}, \link[ASRgenomics]{G.inverse},
#'   \link[ASRgenomics]{kinship.pca}, and \link[ASRgenomics]{snp.pca} that are available
#'   in the package \pkg{ASRgenomics}.
#'
#' @details
#'
#' A summary of the workflow considered in \code{pre.gwas} is:
#' \itemize{
#'   \item 1: Renames markers (if requested).
#'   \item 2: Removes individual records with \code{NA} on \code{resp} and \code{weights}.
#'   \item 3: Organizes variables/columns in \code{pheno.data} (if necessary).
#'   \item 4: Pre-matches \code{pheno.data} and \code{geno.data}.
#'   \item 5: Runs \link[ASRgenomics]{qc.filtering} on \code{geno.data} to obtain filtered marker matrix.
#'   \item 6: Rounds genotypic values on \code{geno.data} (if decimals are present) as required by \link[ASRgenomics]{G.matrix}.
#'   \item 7: Obtains the matrix \eqn{\boldsymbol{K}} using the function \link[ASRgenomics]{G.matrix}, or uses the one provided by the user.
#'   \item 8: Runs \link[ASRgenomics]{G.tuneup} on the \eqn{\boldsymbol{K}} matrix with blending, bending, and aligning (if requested by the user).
#'   \item 9: Performs diagnostics on the \eqn{\boldsymbol{K}} matrix using \link[ASRgenomics]{kinship.diagnostics}.
#'   \item 10: Repeats the matching process on all datasets after the \eqn{\boldsymbol{K}} matrix quality control.
#'   \item 11: Reorders all matrices based on \code{levels(pheno.data[[indiv]])}.
#'   \item 12: Obtains the inverse of the \eqn{\boldsymbol{K}} matrix (\eqn{\boldsymbol{Kinv}}) using \link[ASRgenomics]{G.inverse}.
#'   \item 13: If \eqn{\boldsymbol{Kinv}} is ill-conditioned (changes are cumulative) then:
#'     \itemize{
#'       \item 13.1: Runs \link[ASRgenomics]{G.tuneup} on the \eqn{\boldsymbol{K}} matrix with blending of 0.05 using an identity matrix. If the inverse is still ill-conditioned then go to 13.2; otherwise go to 14.
#'       \item 13.2: Runs \link[ASRgenomics]{G.tuneup} on the \eqn{\boldsymbol{K}} matrix again with blending of 0.10 using an identity matrix. If the inverse is still ill-conditioned then stop; otherwise go to 14.
#'     }
#'   \item 14: Obtains matrix \eqn{\boldsymbol{Q}} based on \code{geno.data} or final \eqn{\boldsymbol{K}} matrix.
#' }
#'
#' The majority of the default values used for arguments are those suggested
#' in \pkg{ASRgenomics}, but these can be modified via the ellipsis argument (\code{...}).
#' Nevertheless, some defaults have been changed for better implementation of
#' GWAS downstream procedures.
#'
#' The modified defaults for GWAS are: \code{ncp = 20}, \code{maf = 0.05},
#' \code{marker.callrate = 0.2}, \code{ind.callrate = 0.2}, \code{clean.diagonal
#' = TRUE}, \code{clean.duplicate = TRUE}, \code{diagonal.thr.large = 1.8}, and
#' \code{diagonal.thr.small = 0.2}. It might be important to adjust these
#' parameters along with GWAS runs to identify best analytical approach.
#'
#' The \eqn{\boldsymbol{Q}} matrix can be calculated based on
#' \eqn{\boldsymbol{K}} or \code{geno.data}. Nevertheless, if missing values are
#' still present in \code{geno.data} after cleaning, \code{Q.method = "K"} will
#' be always used, as this method allows for missing values.
#'
#' For more information on the aforementioned functions, please refer to the
#' \pkg{ASRgenomics} help.
#'
#' @return A list with several data frames and objects to be used in downstream
#'   GWAS model fits.
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
#'
#' @export pre.gwas
#'
#' @examples
#' \dontrun{
#' # Check additional arguments from ASRgenomics (parsed with ...).
#' pre.gwas()
#'
#' # Prepare Apricot datasets (simple run).
#' gwas.data <- pre.gwas(
#'   pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
#'   geno.data = geno.apricot)
#'
#' # Prepare Apricot datasets (adding more filters/procedures).
#' gwas.data <- pre.gwas(
#'   pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
#'   weights = "Fructose", geno.data = geno.apricot, map.data = map.apricot,
#'   Q.method = "geno.data",
#'   # Further arguments.
#'   maf = 0.05, marker.callrate = 0.40, pruning.thr = 0.95, by.chrom = TRUE,
#'   clean.diagonal = FALSE, impute = TRUE)
#' ls(gwas.data)
#' gwas.data$plot.scree
#' head(gwas.data$eigenvalues)
#' }
pre.gwas <- function(
    pheno.data = NULL,
    indiv = NULL,
    geno.data = NULL,
    resp = NULL,
    weights = NULL,
    map.data = NULL,
    rename.markers = TRUE,
    Q.method = "K",
    K = NULL,
    return = c("pheno.data", "geno.data", "K", "Kinv", "map.data", "Q",
               "plot.scree", "eigenvalues"),
    message = TRUE,
    ...
)
{
  if (is.null(indiv) &
      is.null(geno.data) &
      is.null(map.data) &
      is.null(pheno.data)) {
    return(asrgen.args())
  }
  ellipsis <- list(...)
  if (is.null(ellipsis$ncp)){
    ellipsis$ncp = 20
  }
  if (is.null(ellipsis$maf)){
    ellipsis$maf = 0.05
  }
  if (is.null(ellipsis$marker.callrate)){
    ellipsis$marker.callrate = 0.2
  }
  if (is.null(ellipsis$ind.callrate)){
    ellipsis$ind.callrate = 0.2
  }
  if (is.null(ellipsis$clean.diagonal)){
    ellipsis$clean.diagonal = TRUE
  }
  if (is.null(ellipsis$clean.duplicate)){
    ellipsis$clean.duplicate = TRUE
  }
  if (is.null(ellipsis$diagonal.thr.large)){
    ellipsis$diagonal.thr.large = 1.8
  }
  if (is.null(ellipsis$diagonal.thr.small)){
    ellipsis$diagonal.thr.small = 0.2
  }
  if (is.null(ellipsis$pruning.thr)){
    ellipsis$pruning.thr = FALSE
  }
  structure(
    pre.generic(
      pheno.data = pheno.data,
      indiv = indiv,
      resp = resp,
      weights = weights,
      geno.data = geno.data,
      map.data = map.data,
      rename.markers = rename.markers,
      K = K,
      return = return,
      message = message,
      Q.method = Q.method,
      caller = "gwas",
      ellipsis = ellipsis
    )[return],
    class = "gwas.data"
  )
}
