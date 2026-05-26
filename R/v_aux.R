#' Identifies the structure (correlation type) and heterogeneity/homogeneity of variance [adapted from asremlPlus]
#'
#' @param var A character with the name of the factor inside the current term (\emph{e.g.}, \code{family}) (no default).
#' @param term A character with the name of the term (single factor or combination of factors) in
#' \code{random=~} (\emph{e.g.}, \code{family:vm(fishid, Kinv.tmp$fishid)}) (no default).
#' @param G.param The full object named \code{G.param} from the object of class \code{asreml} (no default).
#' @param specials A character vector with the structures allowed (\emph{e.g.}, \code{c("id", "diag")}) (no default).
#' @param residual If \code{TRUE} it means we are working with \eqn{\strong{R}} structures (default = \code{FALSE}).
#'
#' @details
#' It will identify the correlation structure (\emph{e.g.}, \code{"id"}) of the current var within term.
#' Also, will identify if the structure has homogeneous or heterogeneous variance components based
#' on the end of the structure, which can be \code{"v"} or \code{"h"}.
#'
#' @return
#' A list with the structure name and the final (\code{"v"} or \code{"h"}).
#' The output is required by \link{g.asreml} for the reproduction of the structures in matrix form.
#'
#' @author Chris Brien; VSNi
#'
#' @keywords internal
check.special <- function (var, term, G.param, specials, residual = FALSE){
  kspecial <- G.param[[term]][[var]]$model
  if (kspecial == "diag")
    kspecial <- "idh"
  if (kspecial == "dev" | kspecial == "grp")
    kspecial <- "idv"
  nfinal <- nchar(kspecial)
  final <- substr(kspecial, start = nfinal, stop = nfinal)
  if ((final == "v" | final == "h") & kspecial != "sph" &
      kspecial != "dev")
    cortype <- substr(kspecial, start = 1, stop = nfinal - 1)
  else cortype <- kspecial
  if (!(cortype %in% specials)) {
    formula = "random"
    if (residual)
      formula <- "residual"
    stop(paste("No provision has been made for function ",
               kspecial, " in the ", formula, " formula.", sep = ""))
  }
  return(list(cortype = cortype, final = final))
}
#' Gets incidence matrices to estimate \eqn{\strong{G}} [adapted from asremlPlus]
#'
#' @description Collects information from \code{G.param} and produce the correlation
#' matrix used to get \eqn{\strong{G}}.
#'
#' @param var A character with the name of the factor inside the current term (\emph{e.g.}, \code{family}) (no default).
#' @param term A character with the name of the term (single factor or combination of factors) in
#' \code{random=~} (\emph{e.g.}, \code{family:vm(fishid, Kinv.tmp$fishid)}) (no default).
#' @param G.param The full object named \code{G.param} from the object of class \code{asreml} (no default).
#' @param kspecial The output from the \code{check.special} function (no default).
#' @param cond.fac NOT REQUIRED FOR THE CURRENT ACCEPTED CORRELATION STRUCTURES.
#'
#' @return
#' A "correlation matrix" or \eqn{\strong{G}} (if \code{$model} is \code{idv}) associated to a given var (factor).
#'
#' @details
#'
#' If \code{$variance$model == "id"} it means that the variance component is in the
#' \code{$-name of var-$initial} object. In this case the correlation matrix is multiplied by
#' the variance to get \eqn{\strong{G}} already (otherwise it is done outside the function). This
#' happens, for example, in interactions.
#'
#' For now, only \code{id} correlation structures is accepted, others can be added by un-commenting
#' the commented lines (more has to be changed outside).
#'
#' @author Chris Brien; VSNi
#'
#' @keywords internal
g.asreml <- function (var, term, G.param, kspecial, cond.fac = ""){
  if (cond.fac != "")
    cond.fac <- paste(cond.fac, "_", sep = "")
  G <- switch(
    kspecial$cortype,
    id = G.id(
      var = var,
      term = term,
      G.param = G.param,
      cond.fac = cond.fac
    )
  )
  if (kspecial$final == "v") {
    if (kspecial$cortype == "id")
      vpname <- paste(cond.fac, term, "!", var, sep = "")
    else vpname <- paste(cond.fac, term, "!", var, "!var", sep = "")
    if (!is.na(G.param[[term]][[var]]$initial[vpname]))
      G <- G.param[[term]][[var]]$initial[vpname] * G
  }
  return(G)
}
#' Gets a diagonal matrix based on the number of levels of var within term [adapted from asremlPlus]
#'
#' @param var A character with the name of the factor inside the current term (\emph{e.g.}, \code{family}) (no default).
#' @param term A character with the name of the term (single factor or combination of factors) in
#' \code{random=~} (\emph{e.g.}, \code{family:vm(fishid, Kinv.tmp$fishid)}) (no default).
#' @param G.param The full object named \code{G.param} from the object of class \code{asreml} (no default).
#' @param cond.fac Placeholder. Switched off in this function (no default).
#'
#' @return A diagonal matrix.
#'
#' @author Chris Brien; VSNi
#'
#' @keywords internal
G.id <- function (var, term, G.param, cond.fac = ""){
  if (cond.fac == "" | packageVersion("asreml") < "4.2"){
    G <- diag(rep(1, length(G.param[[term]][[var]]$levels)))
    rownames(G) <- colnames(G) <- G.param[[term]][[var]]$levels
  } else{
    G <- diag(rep(1, length(G.param[[term]][[var]]$uniq)))
    rownames(G) <- colnames(G) <- G.param[[term]][[var]]$uniq
  }
  return(G)
}
#' Direct sum of matrices [adapted from asremlPlus]
#'
#' @param matrices List of matrices to be combined via direct sum (no default).
#'
#' @return The direct sum of the provided matrices.
#'
#' @author Chris Brien; VSNi
#'
#' @keywords internal
r.dsum <- function(matrices){
  nr <- lapply(matrices, nrow)
  nc <- lapply(matrices, ncol)
  m <- sum(unlist(nr))
  n <- sum(unlist(nc))
  dsum <- matrix(0, nrow = m, ncol = n)
  r1 <- r2 <- c1 <- c2 <- 0
  for (i in 1:length(matrices)) {
    r1 <- r2 + 1
    c1 <- c2 + 1
    r2 <- r2 + nr[[i]]
    c2 <- c2 + nc[[i]]
    dsum[r1:r2, c1:c2] <- matrices[[i]]
  }
  return(dsum)
}
