#' Calculates repeatability based on the PEV
#'
#' Internal routine that calculates the repeatability based on
#' the Prediction Error Variance (PEV).
#'
#' @param mod An object of class \code{asreml} (default = \code{NULL}).
#' @param component The name of the variance component to be used in the PEV calculations
#' (\emph{e.g.}, "vm(fishid, Kinv$fishid)") (default = \code{NULL}).
#'
#' @return The repeatability of the requested variance component based on PEV.
#'
#' @details
#' Calculations are done using:
#'
#' \deqn{repeatability = 1 - mean(std.error^2) / \sigma^2_c}
#'
#' where: \eqn{std.error} is the average of the standard errors (SE) of the requested factor's
#' solution and \eqn{\sigma^2_c} is its variance component estimate.
#'
#' @keywords internal
repeatability.pev <- function(mod = NULL, component = NULL){
  if (!is.null(mod)){
    if(!is(mod, 'asreml')) {
      stop("Object `mod` should be of class `asreml`.")
    }
  } else {
    stop('No fitted model provided.')
  }
  vparameters <- mod$vparameters
  var.component <- vparameters[names(vparameters) %in% component]
  if (length(var.component) == 0){
    stop("The `component` argument is not found on the provided `asreml` model.")
  }
  solution <- summary(mod, coef = TRUE)$coef.random
  sol.component <- as.data.frame(solution[grep(component, rownames(solution), fixed = TRUE),])
  repeatability <- 1 - mean(sol.component$std.error^2)/(1*var.component)
  if (repeatability < 0 ) {repeatability = 0}
  return(data.frame(Estimate = repeatability))
}
