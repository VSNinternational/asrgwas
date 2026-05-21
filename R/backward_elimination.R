#' Implementation of GWAS model selection by backward elimination
#'
#' @description
#' Performs backward elimination over a set of provided variables based on their
#' \eqn{p-value} when added to a reference model.
#'
#' @param mod A pre-fit \pkg{asreml} model object. Must follow the conditions of the function
#' (default = \code{NULL}).
#' @param data The data used in the process of fitting the \code{asreml} model with additional
#' variables to be tested as fixed effects (default = \code{NULL}).
#' @param try.variables A character vector with the name of variables from \code{data} to be
#' tested in the backward selection procedure (default = \code{NULL}).
#' #' @param pvalue.thr A numeric value with \eqn{p-value} threshold to identify significant markers
#' for the elimination process. Markers will be eliminated until all have a \eqn{p-value} smaller
#' than the value provided (default = \code{0.1})
#' @param maxiter.update A numeric value indicating the number of iterations to be
#' used in \link[asreml]{update.asreml} (default = \code{3}).
#' @param workspace Specifies the workspace to be used by the \code{asreml} function
#' (default = \code{"1Gb"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list containing:
#'
#' \itemize{
#'   \item \code{mod}: The updated base model fitted with the retained variables.
#'   \item \code{kept.variables}: A vector with the names of variables with \eqn{p}-value < \code{pvalue.thr}.
#'   \item \code{aov}: The ANOVA table of the updated model with the kept variables.
#' }
#'
#' @details
#' The elimination technique starts by adding all \code{try.variables} to the base model
#' as fixed effects and calculating their \eqn{p-value}. The variable with the highest \eqn{p-value}
#' is eliminated and the model is re-fitted. The process is repeated until all retained variables
#' have a \eqn{p-value} larger than the provided \code{pvalue.thr}.
#'
#' If \code{vm()} was used to fit the base model, the relationship matrix used \strong{must} be
#' available in \code{.GlobalEnv} under the same name as used in \code{asreml()}.
#'
#' @keywords internal
#'
#' @examples
#' \dontrun{
#' model <- asreml::asreml(yield ~ Variety*Nitrogen, random = ~ Blocks/Wplots,
#'   na.action = list(y='omit', x='omit'), data = oats)
#'
#' bwd.sel <- backward.elimination(
#'   mod = model, data = oats, try.variables = c("Row", "Column"),
#'   pvalue.thr = 0.1, maxiter.update = 10)
#'
#' bwd.sel <- backward.elimination(
#'   mod = model, data = oats,  try.variables = c("Row", "Column"),
#'   pvalue.thr = 0.1, maxiter.update = 10)
#'
#' bwd.sel$aov
#' bwd.sel$kept.variables
#' }
backward.elimination <- function(
    mod = NULL,
    data = NULL,
    try.variables = NULL,
    pvalue.thr = 0.1,
    maxiter.update = 3,
    workspace = "1Gb",
    message = TRUE
){
  .argument_class(.data = data, .class = c("data.frame"))
  .argument_class(.data = mod, .class = "asreml")
  if (pvalue.thr <= 0 || pvalue.thr >= 1){
    stop("The `pvalue.thr` argument should follow the condition 0 < pvalue.thr < 1.")
  }
  if (maxiter.update < 0){
    stop("The `maxiter.update` argument should follow the condition 0 <= maxiter.")
  }
  .argument_class(.data = message, .class = "logical")
  if (!any(try.variables %in% names(data))){
    stop("Some of values of `try.variable` not present in `data`.")
  }
  original.data.name <- as.list(mod$call)$data
  original.data.name <- deparse(substitute(original.data.name))
  mod$call <- fix.call.data_(mod_ = mod, replacement_ = "data")
  if(any(grepl(pattern = "\\-|\\*|\\/|\\%|\\+|\\:", try.variables, fixed = TRUE))){
    stop("Function `asreml` does not accept variable names with mathematical operators (+ - * / % :).")
  }
  variable.loop <- 1:length(try.variables)
  for(holder in variable.loop) {
    try.variables.call <- paste0(try.variables, collapse = "+")
    cur.call <- paste0(
      "cur.mod<-asreml::update.asreml(mod, fixed=.~.+" ,
      try.variables.call, ", workspace ='", workspace,
      "', maxiter=maxiter.update)")
    .silence(eval(expr = parse(text = cur.call)))
    .silence(
      aov <- asreml::wald.asreml(cur.mod, denDF = 'numeric', ssType = 'conditional')$Wald)
    wald.variables <- try.variables[try.variables %in% rownames(aov)]
    marker.pvals <- aov[wald.variables, "Pr"]
    names(marker.pvals) <- wald.variables
    higher.p <- marker.pvals[which.max(marker.pvals)]
    if(higher.p < pvalue.thr) {
      break
    }
    if (message){
      message("Eliminate variable: ", names(higher.p), " [p-value: ", round(higher.p, 4), "].")
    }
    try.variables <- wald.variables[!wald.variables %in% names(higher.p)]
  }
  cur.mod$call <- fix.call.data_(mod_ = cur.mod, replacement_ = original.data.name)
  names(try.variables) <- NULL
  return(list(mod = cur.mod, aov = aov, kept.variables = try.variables))
}
