#' An internal wrapper to calculate VC and PEV heritability/repeatability
#'
#' @param mod An object of class \code{asreml} (default = \code{NULL}).
#' @param herit.numerators The name of the variance component(s) to be used in the numerator of the
#' calculations (\emph{e.g.}, "vm(fishid, Kinv$fishid)") (default = \code{NULL}).
#' @param skip.vc If FALSE, heritability based on variance components
#' will not be calculated (default = \code{FALSE}).
#'
#' @return The repeatabilities based on the variance components and prediction
#' error variances.
#'
#' @keywords internal
repeatability.wrap <- function(mod = NULL, herit.numerators = NULL, skip.vc = FALSE){
  if (skip.vc){
    message("The h2.vc calculation is not available.")
  }
  h2.vc <- list()
  h2.pev <- list()
  for(gen_i in herit.numerators){
    if(!skip.vc) {
      h2.vc[[gen_i]] <-
        repeatability.vc(mod = mod, numerator = gen_i, bound.exclusions = NULL, name = "h2.vc")
    }
    h2.pev[[gen_i]] <- repeatability.pev(mod = mod, component = gen_i)
  }
  if(!skip.vc){
    if (length(herit.numerators) > 1){
      h2.vc[["joint"]] <-
        repeatability.vc(mod = mod, numerator = herit.numerators, bound.exclusions = "B", name = "h2.vc")
    }
    h2.vc <- as.data.frame(rbindlist(h2.vc))
    if (length(herit.numerators) > 1){
      rownames(h2.vc) <- c(herit.numerators, "joint")
    } else {
      rownames(h2.vc) <- herit.numerators
    }
  } else {h2.vc <- NULL}
  h2.pev <- as.data.frame(rbindlist(h2.pev))
  rownames(h2.pev) <- herit.numerators
  return(list(heritability = list("h2.vc" = h2.vc, "h2.pev" = h2.pev)))
}
