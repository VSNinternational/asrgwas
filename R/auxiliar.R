#' Gets the term name for vms (with relationship matrix)
#'
#' @param factor A character vector to be placed in vm(*here*, source$*here*) (default = \code{NULL}).
#' @param source A character vector to be placed in vm(term, *here*$term) (default = \code{"Kinv"}).
#'
#' @return A character vector, \emph{e.g.}, \code{"vm(genA, Kinv$genA)" "vm(genD, Kinv$genD)"}
#'
#' @keywords internal
#'
#' @examples
#' ASRgwas:::get.vm.termvar(c("genA", "genD"), source = "gwas.data$Kinv")
get.vm.termvar <- function(factor = NULL, source = "Kinv"){
  if(is.null(factor)) return(NULL)
  paste("vm(", factor, ", ", paste0(source, "$", factor, ")"), sep = "")
}
#' Moore-Penrose pseudo-inverse
#'
#' @param G A matrix to be inverted (default = \code{NULL}).
#' @param eig.tol Defines relative relevance (\emph{i.e.}, non-zero) of eigenvalues
#' compared to the largest one. It determines which threshold of eigenvalues
#' will be treated as zero (default = \code{1e-06}).
#'
#' @keywords internal
#'
#' @return A Moore Penrose pseudo-inverse of G matrix.
moore.penrose <- function(G = NULL, eig.tol = 1e-06){
  G.dec <- svd(G)
  rel.thr <- max(0, eig.tol * G.dec$d[1])
  keep <- G.dec$d > rel.thr
  if (all(keep)){
    return(
      G.dec$v %*% (1/G.dec$d * t(G.dec$u))
    )
  } else {
    return(
      G.dec$v[, keep, drop = FALSE] %*% ((1/G.dec$d[keep]) * t(G.dec$u[, keep, drop = FALSE]))
    )
  }
}
#' Performs Woodbury's inverted matrix update for row/column elimination
#'
#' @param X An inverted matrix to be updated (default = \code{NULL}).
#' @param na A vector indicating the indices of rows/columns to be dropped (default = \code{NULL}).
#'
#' @return The updated matrix inverse without dropped rows/columns.
#'
#' @keywords internal
#'
woodbury <- function(X, na) {
  A <- X[-na, -na, drop = FALSE]
  B <- X[-na,  na, drop = FALSE]
  C <- X[ na, -na, drop = FALSE]
  D <- X[ na,  na]
  return(A - B %*% solve(D) %*% C)
}
#' Performs Schulz-type inverted matrix update for row/column elimination
#'
#' @param X A matrix to be inverted (default = \code{NULL}).
#' @param Xinv.init An initial guess of the updated inverse (default = \code{NULL}).
#' @param na A vector indicating the indices of rows/columns to be dropped on \code{X}
#' (default = \code{NULL}).
#' @param niter An integer indicating the number of iterations to carry out
#' (theoretically, the more iterations, the closer the approximated inverse is to the true inverse)
#' (default = \code{2})
#'
#' @return The updated matrix inverse with 0 on dropped rows/columns.
#'
#' @keywords internal
schulz <- function(X = NULL, Xinv.init = NULL, na = NULL, niter = 2){
  update.range <- 1:niter
  I2 <- diag(x = 2, nrow = ncol(X))
  Xinv.init[na, ] <- 0
  Xinv.init[, na] <- 0
  for (i in update.range) {
    Xinv.init <- Xinv.init %*% (I2 - X %*% Xinv.init)
  }
  return(Xinv.init)
}
#' Subset ellipsis based on a functions formals.
#'
#' @description
#' Sometimes we want the ellipsis to bear the input for more than one internal function.
#' Therefore we can split it according to the formals of each function.
#'
#' @param ellipsis The ellipsis (\emph{i.e.}, ...) (default = \code{NULL}).
#' @param fun The function to check formals and collect respective input from ... (default = \code{NULL}).
#'
#' @return A list with the ... corresponding to the function.
#'
#' @keywords internal
ellipsis.subset <- function(ellipsis = NULL, fun = NULL){
  vals <- formals(fun)
  args <- formalArgs(fun)
  vals <- lapply(vals, function(x) if (is.call(x)) return(x[[2]]) else x)
  ellipsis.args <- names(ellipsis)
  vars2get <- match(args, ellipsis.args)
  vars2get <- na.exclude(vars2get)
  vals[args %in% ellipsis.args] <- ellipsis[vars2get]
  return(vals)
}
#' Changes data argument in \pkg{asreml} call
#'
#' @description
#' This function will change the call by adjusting the name to be used in
#' the \code{data} argument.
#'
#' @param mod_ An \pkg{asreml} object (default = \code{NULL}).
#'
#' @return The call with the adjusted values in the \code{data} argument.
#'
#' @keywords internal
fix.call.data_ <- function(mod_ = NULL, replacement_ = NULL){
  call.list <- as.list(mod_$call)
  call.list$data <- str2lang(replacement_)
  return(as.call(call.list))
}
#' Change group argument \pkg{asreml} call to match \code{.$mf}
#'
#' @description
#' This function will correct the call by adjusting the variables to be used in
#' the \code{group} argument. The positions might change from the used data and the data
#' reported in \code{.$mf}.
#'
#' @param mod_ An \pkg{asreml} object that used the \code{group} argument (default = \code{NULL}).
#'
#' @return The call with the adjusted values in the \code{group} argument.
#'
#' @keywords internal
fix.call.grp_ <- function(mod_ = NULL){
    mf <- mod_$mf
  grp.name <- names(attributes(mf)$GROUP)
  if (packageVersion("asreml") < "4.2")
    grp.covs <- attributes(mf)$GROUP[[grp.name]]
  else
    grp.covs <- attributes(mf)$model.terms$fixed$Vars[[paste0("grp(", grp.name, ")")]]$Lvls
  grp.covs.pos <- which(names(mf) %in% grp.covs)
  tmp.call <- mod_$call
  tmp.call$group[[grp.name]] <- grp.covs.pos
  return(tmp.call)
}
#' Removes G.param and R.param arguments from \pkg{asreml} object call
#'
#' @param mod_ An \strong{updated} (with \link[asreml]{update.asreml}) \pkg{asreml} object (default = \code{NULL}).
#'
#' @return The call without the \code{G.param} and \code{R.param} arguments.
#'
#' @keywords internal
fix.call.param_ <- function(mod_ = NULL){
  call.list <- as.list(mod_$call)
  call.list <- call.list[!names(call.list) %in% c("G.param", "R.param")]
  return(as.call(call.list))
}
#' Identifies high correlations in matrix of data
#'
#' @param matrix_ An object of class matrix to be checked for high correlations (default = \code{NULL}).
#' @param cor.thr_ A numeric value with the threshold value to generate warnings (default = 0.95).
#'
#' @return Returns a \code{data.frame} with highly correlated combinations.
#'
#' @keywords internal
correlation.report_ <- function(matrix_ = NULL, cor.thr_ = 0.95){
  matrix.cor <- cor(matrix_, use = "pairwise.complete.obs")
  matrix.cor <- 1 * lower.tri(matrix.cor) * matrix.cor
  matrix.cor <- ASRgenomics::full2sparse(matrix.cor)
  matrix.cor.partial <- matrix.cor[abs(matrix.cor[, "Value"]) > cor.thr_ & matrix.cor[, "Value"] < 1, , drop = FALSE]
  matrix.cor.partial[, "Row"] <- attributes(matrix.cor)$rowNames[as.numeric(matrix.cor.partial[, "Row"])]
  matrix.cor.partial[, "Col"] <- attributes(matrix.cor)$colNames[as.numeric(matrix.cor.partial[, "Col"])]
  colnames(matrix.cor.partial)[3] <- "Correlation"
  matrix.cor.partial <- matrix.cor.partial[complete.cases(matrix.cor.partial),]
  if (length(matrix.cor.partial) > 0){
    warning("Correlations that may impair model fitting (r > ", cor.thr_,
            ") identified between at least one pair of columns.")
  }
  matrix.cor.complete <- matrix.cor[abs(matrix.cor[, "Value"]) == 1, , drop = FALSE]
  matrix.cor.complete[, "Row"] <- attributes(matrix.cor)$rowNames[as.numeric(matrix.cor.complete[, "Row"])]
  matrix.cor.complete[, "Col"] <- attributes(matrix.cor)$colNames[as.numeric(matrix.cor.complete[, "Col"])]
  colnames(matrix.cor.complete)[3] <- "Correlation"
  matrix.cor.complete <- matrix.cor.complete[complete.cases(matrix.cor.complete), , drop = FALSE]
  if (length(matrix.cor.complete) > 0){
    warning("Complete correlation identified between at least one pair of columns, check returned object for more information.")
    return(as.data.frame(matrix.cor.complete))
  } else {
    return(cor.markers = NULL)
  }
}
#' Collect object names from call.
#'
#' @param call A model call (default = \code{NULL})
#'
#' @return A character vector with names of objects present in call.
#'
#' @keywords internal
obj.in.call_ <- function(call = NULL){
  operators <- "\\$|\\,|\\(|\\)"
  call.pieces <- unlist(sapply(X = as.character(call), FUN = strsplit, split = operators))
  names(call.pieces) <- NULL
  call.pieces <-trimws(call.pieces)
  call.pieces <- call.pieces[nchar(call.pieces) != 0]
  call.pieces <- call.pieces[sapply(X = call.pieces, FUN = exists)]
  call.pieces <- unique(call.pieces)
  call.pieces <-
    call.pieces[
      !sapply(X = call.pieces,
              FUN = function(x) is.function(try(getFunction(x), silent = TRUE)))]
  if (length(call.pieces) == 0){
    return(NULL)
  } else {
    return(call.pieces)
  }
}
#' @keywords internal
.check.objects.list <- function(.data = NULL,
                                .class = NULL) {
  data.name_ <- deparse(substitute(.data))
  if (!any(class(.data) %in% .class))
    stop(paste0("Objects inside list \'", data.name_, "' should be of class(es) ", paste0(.class, collapse = " or "), "."))
}
