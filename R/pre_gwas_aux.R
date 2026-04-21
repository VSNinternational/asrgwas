#' \pkg{ASRgenomics} function arguments for \link{pre.gwas}
#'
#' @return NULL
#' @keywords internal
asrgen.args <- function(){
  Argument = c("na.string", "marker.callrate", "ind.callrate", "maf", "heterozygosity", "Fis",
               "impute", "method", "criteria", "pruning.thr", "by.chrom", "window.n",
               "overlap.n", "iterations", "seed", "A", "blend", "pblend",
               "bend", "eig.tol", "align", "diagonal.thr.large", "diagonal.thr.small", "duplicate.thr",
               "clean.diagonal", "clean.duplicate", "rcn.thr", "scale", "ncp")
  Description = c(
    "A character string that represents missing values (default = \"NA\").",
    "Marker call rate threshold (default = 0.2).",
    "Individual call rate threshold (default = 0.2).",
    "MAF (Minor Allele Frequency) threshold (default = 0.05).",
    "Heterozygosity threshold (default = 1, i.e. no removing).",
    "Fis (inbreeding coefficient) threshold (default = 1, i.e. no removing).",
    "If TRUE, mean imputation of missing values is performed (default = FALSE).",
    "Method of GRM (Genomic Relationship Matrix) calculation (\"VanRaden\", \"Yang\", \"Su\", \"Vitezica\") (default = \"VanRaden\").",
    "Criteria (\"callrate\" or \"maf\") to choose marker to drop in pruning (default = \"callrate\").",
    "Correlation threshold to identify redundant markers in pruning (default = NULL; i.e. no pruning).",
    "If TRUE, the pruning is performed independently by chromosome (default = FALSE).",
    "Number of markers to consider in each window to perform pruning (default = 50).",
    "Number of markers to overlap between consecutive windows for pruning (default = 5).",
    "Number of sequential times the pruning procedure should be executed (default = 10).",
    "A seed for reproducibility in pruning (default = NULL).",
    "A pedigree relationship matrix for blending/alignment (default = NULL).",
    "If TRUE, a blending of K is performed (default = FALSE).",
    "The proportion of the pedigree or identity matrix to use in blending (default = NULL).",
    "If TRUE, a bending is performed (default = FALSE).",
    "Determines which threshold of eigenvalues will be treated as zero (default = 1e-06).",
    "If TRUE, the GRM (Genomic Relationship Matrix) is aligned to the matching pedigree matrix (default = FALSE).",
    "A threshold value to flag large diagonal values (default = 1.8).",
    "A threshold value to flag small diagonal values (default = 0.2).",
    "A threshold value to flag possible duplicates (default = 0.95).",
    "If TRUE, filters out values > diagonal.thr.large and < diagonal.thr.small (default = TRUE).",
    "If TRUE, return a kinship matrix without duplicate individuals (default = TRUE).",
    "A threshold for identifying the K matrix as an ill-conditioned matrix (default = 1e-12).",
    "If TRUE, the PCA (Principal Components Analysis) will scale the kinship matrix (default = TRUE).",
    "The number of PC (Principal Components) dimensions to provide in the output data frame (default = 20)."
  )
  cat("\033[34mAdditional ASRgenomics arguments that can be passed to pre.gwas():\033[39m\n")
  print.data.frame(
    data.frame(
      Argument = Argument, Description = Description
    ), right = F
  )
}
