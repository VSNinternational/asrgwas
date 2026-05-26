#' Simple argument-oriented function to build \pkg{asreml} calls
#'
#' @param resp A character with name to be used as response variable in \pkg{asreml} call (default = \code{NULL}).
#' @param fixedf A character vector with variable names to be used as fixed effects in \pkg{asreml} call.
#' Interactions/nesting are allowed (default = \code{NULL}).
#' @param group.cols A numeric vector with the number of the column(s) in data to be grouped with the
#' \code{grp()} function (default = \code{NULL}).
#' @param cov A character vector with variable names to be used as covariates in \pkg{asreml} call
#' (default = \code{NULL}).
#' @param randomf A character vector with variable names to be used as random effects in \pkg{asreml} call.
#' Interactions/nesting are allowed (including with randomf.vm) (default = \code{NULL}).
#' @param randomf.vm A character vector with variable names to be used as random effects with
#' user defined relationship matrix in \pkg{asreml} call.
#' Interactions/nesting NOT allowed. Specify all interactions in \code{randomf} (default = \code{NULL}).
#' @param vm.matrix.name The name of the object containing the relationship matrix to be passed to
#' \code{vm()} (default = \code{NULL}).
#' @param data.name A character to be used as name of the data frame in the \pkg{asreml} call (default = \code{NULL}).
#' @param residual A character with the name of a variable in \code{pheno.data} to be used as
#' conditional factor for specifying heterogeneous residuals by level of this factor.
#' @param weights A character with name of the numeric variable in \code{pheno.data} with the
#' weights of each prediction (\emph{i.e.}, mean estimate) (default = \code{NULL}).
#' Note that \code{pheno.data} will be internally ordered by this factor (default = \code{NULL}).
#' @param family A character indicating the family of distribution for \code{resp}.
#' Options are: \code{"gaussian"} and \code{"binomial"} (default = \code{"gaussian"}).
#'
#' @return A character string with the \pkg{asreml} call.
#'
#' @keywords internal
#'
#' @details
#' Currently available models in this function can be built with:
#' any combination of fixed and random effects;
#' covariates; groups (\code{grp()}); user defined structures; \code{vm():ind()} interactions/nesting terms;
#' heterogeneous residuals with \code{dsum(|)}; addition of weights; Gaussian and binomial families.
#'
#' @examples
#'\dontrun{
#' # A dummy example.
#' asreml.call <-
#'   ASRgwas:::call.generator(asreml.obj.name = "mod",
#'     resp = "resp", fixedf = c("year", "block"), group.cols = c(10,12),
#'     cov = "stand", randomf = c("gen:year"), randomf.vm = "gen",
#'     vm.matrix.name = "Kinverse", data.name = "phenotypic.data",
#'     residual = "year", weights = NULL, family = "gaussian", workspace = "4Gb")
#' str2lang(asreml.call)
#'
#' # An (adapted) example from asreml-r manual.
#' asreml.call <-
#'   ASRgwas:::call.generator(
#'     asreml.obj.name = "mod",
#'     resp = "yield", fixedf = c("Blocks", "Nitrogen"),
#'     cov = "Column", randomf = c("Variety", "Subplots"),
#'     data.name = "oats", family = "gaussian", workspace = "1Gb")
#' eval(parse(text = asreml.call))
#' }
call.generator <- function(
    resp = NULL, fixedf = NULL,
    group.cols = NULL, cov = NULL,
    randomf = NULL, randomf.vm = NULL,
    vm.matrix.name = NULL, data.name = NULL,
    residual = NULL, weights = NULL,
    family = NULL, dispersion = 1, total = NULL,
    workspace = NULL){
  code.asr <- as.character()
  code.asr["fixed"] <- paste0("asreml::asreml(fixed=", resp, "~1")
  code.asr["random"] <- "random=~"
  if (family == "gaussian"){
    if (!is.null(residual)){
      code.asr["residual"] <- paste0("residual=~dsum(~units|", residual, ")")
    }
    if(!is.null(weights)){
      code.asr["residual"] <- "family=asreml::asr_gaussian(dispersion = 1)"
    }
    if(is.null(residual) & is.null(weights)){
      code.asr["residual"] <- "residual=~idv(units)"
    }
  }
  if (family == "binomial"){
    if (is.null(total)){
      code.asr["residual"] <-
        paste0("residual=~id(units),family=asreml::asr_binomial(link=\'logit\',dispersion=",
               dispersion, ")")
    }
    if (!is.null(total)){
      code.asr["residual"] <-
        paste0("residual=~id(units),family=asreml::asr_binomial(link=\'logit\',dispersion=",
               dispersion, ",total=\'", total, "\')")
    }
  }
  code.asr["workspace"] <- paste0("workspace=\'", workspace, "'")
  code.asr["na.action"] <- "na.action=list(x='include',y='include')"
  if (!is.null(weights)){
    code.asr["weight"] <- paste0("weights=", weights)
  }
  if (!is.null(group.cols)){
    code.asr["fixed"] <- paste(code.asr["fixed"], "grp(Q)", sep = "+")
    code.asr["group"] <- paste('group=list(Q=', group.cols[1] , ':', group.cols[2], ')', sep='')
  }
  code.asr["data"] <- paste0("data=", data.name, ")")
  if (!is.null(cov))
    code.asr["fixed"] <- paste(code.asr["fixed"], paste(cov, collapse = "+"), sep = "+")
  if (!is.null(fixedf))
    code.asr["fixed"] <- paste(code.asr["fixed"], paste(fixedf, collapse = "+"), sep = "+")
  if (!is.null(randomf)){
    tmp.rdm.int <- strsplit(x = randomf, split = ":", fixed = TRUE)
    tmp.rdm.inc <- lapply(tmp.rdm.int, function(t) any(t %in% randomf.vm))
    for (t in seq_along(tmp.rdm.inc)) {
      if (isTRUE(tmp.rdm.inc[[t]])){
        vm.inc <- tmp.rdm.int[[t]] %in% randomf.vm
        first.idv.pos <- min(which(!vm.inc == TRUE))
        new.vm.var <- paste("vm(", tmp.rdm.int[[t]][vm.inc], ", ",
                            paste0(vm.matrix.name, "$", tmp.rdm.int[[t]][vm.inc] , ")"), sep = "")
        new.idv.var <- paste0("idv(", tmp.rdm.int[[t]][first.idv.pos], ")")
        tmp.rdm.int[[t]][first.idv.pos] <- new.idv.var
        tmp.rdm.int[[t]][vm.inc] <- new.vm.var
        randomf[t] <- paste0(tmp.rdm.int[[t]], collapse = ":")
      }
    }
    inter.randomf <- paste(randomf, collapse = "+")
  } else {
    inter.randomf <- NULL}
  if(!is.null(randomf.vm)){
    vm.vars <- paste("vm(", randomf.vm, ", ", paste0(vm.matrix.name, "$", randomf.vm , ")"), sep = "")
  } else vm.vars <- NULL
  if (!is.null(randomf.vm) & !is.null(randomf))
    code.asr["random"] <- paste(code.asr["random"], paste(vm.vars, collapse = "+"), inter.randomf, sep = "+")
  else if (!is.null(randomf.vm))
    code.asr["random"] <- paste(code.asr["random"], paste(vm.vars, collapse = "+"), sep = "+")
  else if (!is.null(randomf))
    code.asr["random"] <- paste(code.asr["random"], paste(randomf, collapse = "+"), sep = "+")
  if (is.null(randomf.vm) & is.null(randomf)){
    code.asr <- code.asr[!names(code.asr)=="random"]
  }
  code.asr["fixed"] <- paste(code.asr["fixed"], sep='')
  code.asr <- paste0(code.asr, collapse=',')
  return(code.asr)
}
