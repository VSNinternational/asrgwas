#' Calculates repeatability based on variance components
#'
#' @description
#' Uses \link[asreml]{vpredict} to calculate the repeatability
#' based on the variance components for a given (set of) component(s).
#'
#' @param mod An object of class \code{asreml} (default = \code{NULL}).
#' @param numerator The name of the variance component(s) to be used in the numerator of the
#' calculations (\emph{e.g.}, "vm(fishid, Kinv$fishid)") (default = \code{NULL}).
#' @param bound.exclusions The code of the boundary to be excluded (\emph{e.g.}, c("F", "B");
#' default = \code{NULL})
#'
#' @return A data frame with the estimates, standard error, and formula for the h2.vc calculation.
#'
#' @details
#' The algorithm uses the values provided in \code{numerator} as the numerator in the formula.
#' All variance components are added to the denominator, except when its bound conditions match values
#' passed to \code{bound.exclusions}.
#'
#' The formula can be represented as:
#'
#' \deqn{repeatability = \sigma^2_n / \sigma^2_d}
#'
#' where: \eqn{\sigma^2_n} is the sum of the variance components passed to the \code{numerator} argument,
#' and \eqn{\sigma^2_d} is the sum of all variance components.
#'
#' If there are heterogeneous residual variances, these will be averaged before calculations.
#'
#' @keywords internal
repeatability.vc <- function(
  mod = NULL,
  numerator = NULL,
  bound.exclusions = NULL,
  name = "repeatability"
) {
  if (!is.null(mod)) {
    if (!is(mod, 'asreml')) {
      stop("The object `mod` should be of class `asreml`.")
    }
  } else {
    stop('No fitted model provided.')
  }
  vparameters <- mod$vparameters
  bounded.vars <- which(summary(mod)$varcomp$bound %in% bound.exclusions)
  cor.vars <- grep(pattern = "cor", x = names(vparameters))
  units.excl.R.var <- grep(pattern = "units!R", x = names(vparameters))
  excl.R.vars <- grep(pattern = "!R", x = names(vparameters))
  if (length(excl.R.vars) > 1) {
    residual.vars <- excl.R.vars
  } else {
    residual.vars <- NULL
  }
  numerator.vars <- c()
  for (g in seq_along(numerator)) {
    numerator.vars <- append(
      numerator.vars,
      tryCatch(
        expr = grep(
          pattern = numerator[g],
          x = names(vparameters),
          fixed = TRUE
        ),
        error = function(holder) {
          return(NULL)
        }
      )
    )
  }
  numerator.vars <- unique(numerator.vars)
  numerator.inter.vars <- grep(":", names(vparameters))
  if (length(numerator.inter.vars) > 0) {
    numerator.vars <- na.exclude(numerator.vars[
      numerator.vars != numerator.inter.vars
    ])
  }
  other.vars <- seq_along(vparameters)[na.exclude(
    -c(bounded.vars, residual.vars, cor.vars, numerator.vars, units.excl.R.var)
  )]
  if (length(other.vars) > 0) {
    other.vars.formula <- paste0("V", other.vars, collapse = "+")
  } else {
    other.vars.formula <- NULL
  }
  numerator.formula <- paste0("V", numerator.vars, collapse = "+")
  if (length(residual.vars) > 1) {
    het.res.formula <- paste0(
      "(",
      paste0("V", residual.vars, collapse = "+"),
      ")/",
      length(residual.vars)
    )
  } else {
    het.res.formula <- NULL
  }
  phen.formula <- paste0(
    c(numerator.formula, other.vars.formula, het.res.formula),
    collapse = "+"
  )
  herit.call <- paste0(name, "~(", numerator.formula, ")/(", phen.formula, ")")
  herit.val <- asreml::vpredict(object = mod, xform = as.formula(herit.call))
  return(data.frame(herit.val, Formula = herit.call))
}
