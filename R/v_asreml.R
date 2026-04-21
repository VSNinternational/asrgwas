#' Calculates the variance-covariance of \eqn{y} (\eqn{\strong{V}} matrix) [adapted from asremlPlus]
#'
#' @description
#' This function calculates variance-covariance of \eqn{y} based on an object of class \code{asreml}
#' (fitted model).
#'
#' We are not using asremlPlus here because we are allowing interactions with vms(). It would be
#' difficult to get the additional matrices outside the function. Lots of things have been commented
#' out but kept here for future reference and comparison to the main function.
#'
#' @param asreml.obj An \pkg{asreml} object (default = \code{NULL}).
#' @param which.matrix This will probably be removed.
#' @param Ks A list with (non-inverted) full relationship matrices. See details (default = \code{NULL}).
#' @param weights A character indicating a variable in \code{asreml.obj} containing
#' residual weights (default = \code{NULL}).
#'
#' @details
#'
#' Changes from the original function (asremlPlus):
#'
#' * \code{vm} structure has been added;
#' * interactions/nesting with \code{vm} has been added;
#' * identification of \code{gammaPar} has changed;
#' * added possibility to include residual weights in the calculations;
#' * several things have been switched off for now;
#'
#' Terminology:
#'
#' \code{term}: a composite model term; e.g. \code{vm(ID, Kinv):Env}
#'
#' A \code{ranterm} is anything passed to \code{random=~} in \pkg{asreml} (e.g. \code{genotype:env}).
#' This is later split into vars (e.g. if term is \code{genoype:env}, we will get
#' \code{genotype} and \code{env} as vars).
#'
#' \code{termvars}: are the vars in each term. Might be one or more.
#' Argument \code{Ks} should be passed as follows: assuming \eqn{\strong{KA}} and \eqn{\strong{KD}} are two
#' relationship matrices inverted to \eqn{\strong{KinvA}} and \eqn{\strong{KinvD}}, respectively, which were used as
#' \code{vm} structures for fitting \code{asreml.obj}; the \code{Ks} should be
#' \code{list("vm(genA, KinvA)" = KA, "vm(genD, KinvD)" = KD)}, where \code{vm(genA, KinvA)} and
#' \code{vm(genD, KinvD)} are the names of the random terms; and \code{KA} and \code{KD}
#' are the original relationship matrices.
#'
#' Additional information:
#'
#' * Structures accepted for the random terms are: \code{id}, \code{idv}, and \code{vm}.
#' More can be added later on.
#' * Structures accepted for the residual terms are: \code{id}, \code{idv}, and \code{dsum}.
#' More can be added later on.
#' * Interaction/nesting with two \code{vm} is not permitted by \pkg{asreml} (too complex).
#' * It iterates through each term's variable (factor).
#' At each run, the Kronecker will be calculated with the previous terms.
#' This is used when there are interactions (more than one variable in the term).
#' * G.param has information about each term. Inside each term we have a \code{$varinace}
#' object and \code{$-name of var-} objects. The variance of the term might be on
#' \code{$varinace$initial} or in one of the vars \code{$-name of var-$initial}.
#' We also have to collect the \code{$model} associated to the object. If \code{model == "idv"};
#' the \code{$initial} will have the variance component. If not \code{"idv"}
#' it will be a correlation matrix, e.g. \code{vm}, \code{id}, etc.
#' * When \code{vm} is used, the relationship matrix is the correlation matrix itself.
#' It can be directly multiplied by the associated variance to get \eqn{\strong{G}}.
#' * There is a loop to sum up the variance of all random terms, and get \eqn{\strong{G}}.
#' \eqn{\strong{G}} (variance-covariance of random terms) and \eqn{\strong{R}}
#' (variance-covariance of residuals) are estimated separately and then summed up.
#' * If \code{gammaPar} has been used, the \eqn{\strong{V}} matrix is multiplied by the
#' scaling factor \code{$sigma2} at the end.
#'
#' @return The variance-covariance of \eqn{y}: \eqn{\strong{V}} matrix.
#'
#' @author Chris Brien; VSNi
#'
#' @keywords internal
v.asreml <- function (
    asreml.obj = NULL,
    which.matrix = c("V", "G", "R"),
    Ks = NULL,
    weights = NULL)
{
  which.matrix <- match.arg(which.matrix)
  ran.specials <- c()
  res.specials <- NULL
  common.specials <- c("id")
  call <- asreml.obj$call
  data <- asreml.obj$mf
  n.samples <- nrow(data)
  V <- matrix(0, nrow = n.samples, ncol = n.samples)
  G.param <- asreml.obj$G.param
  ranterms <- names(G.param)
  R.param <- asreml.obj$R.param
  resterms <- names(R.param)
  if (which.matrix %in% c("V", "G")) {
    for (term in ranterms) {
      G.param.vars <- G.param[[term]][2:length(G.param[[term]])]
      idv.only.vars <- lapply(G.param.vars, function(x) x$model == "idv")
      if (all(unlist(idv.only.vars)) & length(idv.only.vars) > 1)
        stop("The term `", term, "` has competing variance components.")
      id.only.vars <- lapply(G.param.vars, function(x) x$model == "id")
      id.only.vars <- all(unlist(id.only.vars))
      variance.idv <- G.param[[term]]$variance$model == "idv"
      rm(G.param.vars)
      if (
        variance.idv & id.only.vars) {
        Z <- model.matrix(as.formula(paste("~ - 1 + ", term)), data = data)
        V <- V + G.param[[term]]$variance$initial * Z %*% t(Z)
        rm(id.only.vars, variance.idv, idv.only.vars)
      }
      else {
        termvars <- names(G.param[[term]])[-1]
        cond.fac <- ""
        if (G.param[[term]]$variance$model == "idv")
          G <- G.param[[term]]$variance$initial
        else G <- 1
        term.has.vm <- length(grep(pattern = "vm", x = term)) > 0
        if (term.has.vm){
          clean.term <- lapply(G.param[[term]][2:length(G.param[[term]])], function(x) x$model)
          clean.term <- lapply(clean.term, function(x) names(x))
          clean.term <- unlist(clean.term)
          clean.term <- paste0(clean.term, collapse = ":")
        }
        for (var in termvars) {
          var.has.vm <- length(grep(pattern = "vm", x = var)) > 0
          if (var.has.vm){
            vm.var <- lapply(G.param[[term]][2:length(G.param[[term]])], function(x) x$model)
            vm.var <- vm.var[[which(unlist(lapply(vm.var, function(x) x == "vm")))]]
            vm.var <- names(vm.var)
          }
          if (!var.has.vm){
            kspecial <- check.special(var = var, term = term,
                                      G.param = G.param,
                                      specials = c(ran.specials,
                                                   common.specials),
                                      residual = FALSE)
            G <- kronecker(G, g.asreml(var = var, term = term,
                                       G.param = G.param, kspecial = kspecial,
                                       cond.fac = cond.fac))
          }
          if (var.has.vm){
            if (!identical(levels(data[[vm.var]]), rownames(Ks[[vm.var]])))
              stop("The order of the provided var-cov matrix does not match levels of the variable used in vm().")
            G <- kronecker(G, Ks[[vm.var]])
          }
          if (!exists("tmp.names"))
            if (var.has.vm)
              tmp.names <- paste0(vm.var, levels(data[[vm.var]]))
          else  tmp.names <- paste0(var, G.param[[term]][[var]]$levels)
          else if (var.has.vm)
            tmp.names <- paste0(
              rep(x = tmp.names, each = nlevels(data[[vm.var]])),
              ":", vm.var, levels(data[[vm.var]]))
          else tmp.names <- paste0(
            rep(x = tmp.names, each = length(G.param[[term]][[var]]$levels)),
            ":", var, G.param[[term]][[var]]$levels)
        }
        if (term.has.vm)
          tmp.term <- clean.term
        else  tmp.term <- term
        Z <- model.matrix(as.formula(paste("~ -1 +", tmp.term)), data = data)
        rownames(G) <- colnames(G) <- tmp.names ; rm(tmp.names)
        G <- G[colnames(Z), colnames(Z)]
        V <- V + Z %*% G %*% t(Z)
      }
    }
  }
  if (which.matrix %in% c("V", "R")) {
    nosections <- length(R.param)
    Rlist <- vector(mode = "list", length = nosections)
    names(Rlist) <- resterms
    if (nosections > 1) {
      cond.fac <- strsplit(labels(asreml.obj$formulae$residual)[1],
                           split = "|", fixed = TRUE)[[1]][2]
      if (grepl(",", cond.fac, fixed = TRUE))
        cond.fac <- strsplit(cond.fac, ",", fixed = TRUE)[[1]][1]
      if (grepl(")", cond.fac, fixed = TRUE))
        cond.fac <- strsplit(cond.fac, ")", fixed = TRUE)[[1]][1]
      cond.fac <- trimws(cond.fac)
    }
    else {
      cond.fac <- ""
    }
    for (term in resterms) {
      termvars <- names(R.param[[term]])[-1]
      if (attributes(data)$scale.fit == "sigma")
        R <- R.param[[term]]$variance$initial
      else R <- 1
      for (var in termvars) {
        kspecial <-
          check.special(
            var = var, term = term,
            G.param = R.param,
            specials = c(res.specials, common.specials),
            residual = TRUE)
        R <-
          kronecker(
            R, g.asreml(
              var = var,
              term = term,
              G.param = R.param,
              kspecial = kspecial,
              cond.fac = cond.fac))
      }
      Rlist[[term]] <- R
    }
    if (is.null(weights)) {
      if (length(Rlist) == 1) {
        V <- V + Rlist[[1]]
      }
      else V <- V + r.dsum(Rlist)
    } else {
      V <- V + diag(1/data[[weights]])
    }
  }
  if (attributes(data)$scale.fit == "gamma")
    V <- asreml.obj$sigma2 * V
  return(V)
}
