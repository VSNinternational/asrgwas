#' Implements GWAS model fit using \link[asreml]{update.asreml}
#'
#' This functions runs conventional or P3D GWAS using \link[asreml]{update.asreml}.
#' One marker will be added to the original model at a time and the model updated.
#'
#' @param mod A pre-fit \pkg{asreml} model object. Must follow the conditions of the function
#' (default = \code{NULL}).
#' @param geno.data A matrix with SNP data of form \eqn{n \times p},
#' with \eqn{n} individuals and \eqn{p} markers. Individual and marker names are
#' assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as: 0, 1, 2 (default = \code{NULL}).
#' @param maxiter.update A numeric value indicating the number of iterations to be
#' used in \link[asreml]{update.asreml} under the non-P3D method (default = \code{1}).
#' @param P3D If \code{TRUE} the "Population Parameters Previously Determined" algorithm will be used.
#' If \code{FALSE} the variance components are allowed to vary while fitting each marker individually.
#' (default = \code{FALSE}).
#' @param threads An integer with the number of threads to be used in parallel processing
#' (default = \code{1}).
#' @param obj.parallel A character vector indicating which objects passed to the function that must
#' be exported to clusters (\emph{e.g.}, objects passed to \code{Kinv} argument). If nothing is passed, the
#' algorithm will try to identify the required object. This is only required in Windows as Unix systems
#' share environments by default (default = \code{NULL}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with \code{effect}, \code{z.ratio} and \code{p.value} of each marker present
#' in \emph{geno.data}.
#'
#' @keywords internal
mc.marker.update <- function(mod = NULL, geno.data = NULL, maxiter.update = NULL,
                             threads = 1, obj.parallel = NULL, P3D = FALSE, message = TRUE) {
  asreml::asreml.options(design = FALSE, trace = FALSE)
  family <- attributes(mod$mf)$traits$family
  pheno.data <- mod$mf
  mod$call <- fix.call.data_(mod_ = mod, replacement_ = "pheno.data")
  if (!is.null(mod$call$group)){
    mod$call <- fix.call.grp_(mod_ = mod)
  }
  n.markers <- ncol(geno.data)
  if (P3D){
    if (message) {
      message("Fixing variance components in the base model (P3D).")
    }
    mod$G.param <-
      lapply(mod$G.param, function(fl)
        lapply(fl, function(sl){
          sl$con <- "F" ; return(sl)}
        )
      )
    mod$R.param <-
      lapply(mod$R.param, function(fl)
        lapply(fl, function(sl){
          sl$con <- "F" ; return(sl)}
        )
      )
  }
  if (threads > detectCores()){
    if (message){
      message(
        paste0("The `threads` is larger than available threads, reseting to maximum (",
               detectCores(),  ")."))
      threads <- detectCores()
    }
  }
  if (threads > 1){
    if (message){
      message("Parallel processing selected. Progress assessment only if `threads = 1`.")
    }
  }
  if (threads == 1 & message){
    cli_progress_bar(name = "Marker fitting", total = n.markers, clear = FALSE,
                     format = "Progress: {cli::pb_bar} [{cli::pb_current}/{cli::pb_total}] [{round(cli::pb_rate_raw, 0)} markers/s] [{cli::pb_eta_str}] [elapsed: {cli::pb_elapsed}]")
  }
  iter <- 1:n.markers
  if (message){
    message("Iterating through markers...")
  }
  marker <- NULL
  marker.update_ <- function(marker.index_ = NULL, mod_ = NULL, pheno.data = NULL,
                             maxiter.update_ = NULL, geno.data_ = NULL) {
    marker <<- geno.data_[, marker.index_]
    .silence(
      mod.updated <-
        asreml::update.asreml(
          object = mod_,
          fixed = . ~ . + marker,
          maxiter = maxiter.update_,
          na.action = list(y = 'include', x = 'include'))
      )
    solution <-
      summary(mod.updated, coef = TRUE)$coef.fixed["marker",]
    return(
      list(effect = solution[1],
           std.error = solution[2],
           z.ratio = solution[3],
           ndf = mod.updated$nedf))
  }
  if (threads == 1){
    marker.stats <- list()
    time.start <- proc.time()["elapsed"]
    for (mi in iter){
      marker.stats[[mi]] <- marker.update_(marker.index_ = mi, mod_ = mod, pheno.data = pheno.data,
                                           maxiter.update_ = maxiter.update, geno.data_ = geno.data)
      if (message) cli_progress_update()
    }
  }
  if (threads > 1){
    if (Sys.info()["sysname"] == "Windows"){
      cl <- makeCluster(threads)
      if (is.null(obj.parallel)){
        obj.parallel <- obj.in.call_(call = mod$call)
      }
      clusterExport(cl = cl, obj.parallel)
      time.start <- proc.time()["elapsed"]
      marker.stats <- parLapply(cl = cl, iter, marker.update_, mod_ = mod, pheno.data = pheno.data,
                                          maxiter.update_ = maxiter.update, geno.data_ = geno.data)
      stopCluster(cl)
    }
    if (Sys.info()["sysname"] != "Windows"){
      time.start <- proc.time()["elapsed"]
      marker.stats <- mclapply(iter, marker.update_, mod_ = mod, pheno.data = pheno.data,
                                         maxiter.update_ = maxiter.update, geno.data_ = geno.data,
                                         mc.cores = threads)
    }
  }
  time.finish <- proc.time()["elapsed"] - time.start
  marker.stats <- rbindlist(marker.stats)
  if (family == "gaussian"){
    p.value <- exp(
      log(2) +
        pt(abs(marker.stats$z.ratio), lower.tail = FALSE, log.p = TRUE, df = marker.stats$ndf)
    )
  }
  if (family == "binomial"){
    p.value <- exp(
      log(2) +
        pnorm(abs(marker.stats$z.ratio), mean = 0, sd = 1, lower.tail = FALSE, log.p = TRUE)
    )
  }
  return(list(effect = marker.stats$effect, std.error = marker.stats$std.error,
              z.ratio = marker.stats$z.ratio, p.value = p.value,
              time.taken = time.finish
              ))
}
