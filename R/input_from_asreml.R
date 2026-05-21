#' Collects data from an asreml object
#'
#' @param mod An object of class \code{asreml} (default = \code{NULL}).
#' @param parent.frame The environment in which the data will be assigned (default = \code{1} which
#' means one environment above).
#'
#' @return The \code{pheno.data}, \code{resp}, \code{family}, \code{weight}, \code{residual}
#' on the selected environment.
#'
#' @keywords internal
input.from.asreml <- function(mod = NULL, parent.frame = 1){
  pheno.data <- mod$mf
  attrb <- attributes(pheno.data)
  resp <- attrb$traits$lhs
  family <- attrb$traits$family
  weight <- attrb$weights
  if (length(mod$R.param) > 1){
    labels(mod$formulae$residual)
    residual <- strsplit(labels(mod$formulae$residual)[1],
                         split = "|", fixed = TRUE)[[1]][2]
    if (grepl(")", residual, fixed = TRUE)){
      residual <- strsplit(residual, ")", fixed = TRUE)[[1]][1]
    }
    residual <- trimws(residual)
  } else {
    residual <- NULL
  }
  assign(x = "pheno.data", value = pheno.data, envir = parent.frame(parent.frame))
  assign(x = "resp", value = resp, envir = parent.frame(parent.frame))
  assign(x = "family", value = family, envir = parent.frame(parent.frame))
  assign(x = "weight", value = weight, envir = parent.frame(parent.frame))
  assign(x = "residual", value = residual, envir = parent.frame(parent.frame))
}
