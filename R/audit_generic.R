#' Checks/audits datasets and matrices prior to GBLUP/GWAS
#'
#' @description
#' See \link[ASRgwas]{audit.gwas}
#' @keywords internal
audit.generic <- function(
  pheno.data = NULL,
  indiv = NULL,
  resp = NULL,
  Kinv = NULL,
  geno.data = NULL,
  map.data = NULL,
  Q = NULL,
  caller = c("gwas", "gblup"),
  message = TRUE
) {
  out <- list()
  .argument_class(.data = pheno.data, .class = c("data.frame"))
  .variable_class(
    .data = pheno.data,
    .mandatory = TRUE,
    .variable = indiv,
    .class = "factor",
    .class.action = "message",
    .mutate = TRUE,
    .message = message
  )
  .variable_class(
    .data = pheno.data,
    .mandatory = FALSE,
    .variable = resp,
    .class = "numeric",
    .class.action = "message",
    .mutate = TRUE,
    .message = message
  )
  .argument_class(.data = geno.data, .class = c("matrix", "NULL"))
  .argument_class(.data = Kinv, .class = c("list", "matrix", "NULL"))
  if (is.matrix(Kinv)) {
    Kinv <- list(Kinv)
  }
  ifelse(
    test = is.null(geno.data),
    yes = {
      geno.data.type <- "NIL"
    },
    no = {
      geno.data.type <- "OBJ"
    }
  )
  match.data(
    pheno.data = pheno.data,
    indiv = indiv,
    geno.data = switch(geno.data.type, NIL = NULL, OBJ = list(geno.data)),
    Kinv = Kinv,
    Q = Q,
    map.data = map.data,
    caller = caller,
    message = message
  )
  if (!is.null(geno.data)) {
    miss.geno.data <- sum(is.na(geno.data)) / length(geno.data) * 100
    if (message & miss.geno.data > 0) {
      message(
        "Missing data found in `geno.data` (",
        round(miss.geno.data, 2),
        "%)."
      )
    }
    rm(miss.geno.data)
  }
  if (!is.null(resp)) {
    out$gen.stats <- .by.stats(.data = pheno.data, .index = indiv, .var = resp)
    out$trial.stats <- .by.stats(.data = pheno.data, .index = NULL, .var = resp)
    out$trial.stats <- out$trial.stats[, -1]
  }
  if (!is.null(geno.data)) {
    out$n.genotypes <- nrow(geno.data)
  } else if (!is.null(Kinv)) {
    out$n.genotypes <- length(attributes(Kinv[[1]])$rowNames)
  } else if (!is.null(pheno.data)) {
    out$n.genotypes <- nlevels(pheno.data[[indiv]])
  } else {
    out$n.genotypes <- "unknown"
  }
  if (message) {
    message(paste0(
      "A total of ",
      out$n.genotypes,
      " unique genotypes were found."
    ))
  }
  if (!is.null(geno.data)) {
    out$n.markers <- ncol(geno.data)
    if (message) {
      message(paste0("A total of ", out$n.markers, " markers were found."))
    }
  }
  if (!is.null(geno.data)) {
    out$ind.callrate.stats <- ASRgenomics:::callrate(
      M = geno.data,
      margin = "row"
    )
    out$ind.callrate.stats <- round(range(out$ind.callrate.stats), 2)
    if (message) {
      message(paste0(
        "Genotype call rate range in `geno.data` is: ",
        paste0(out$ind.callrate.stats, collapse = " ~ "),
        "."
      ))
    }
    out$marker.callrate.stats <- ASRgenomics:::callrate(
      M = geno.data,
      margin = "col"
    )
    out$marker.callrate.stats <- round(range(out$marker.callrate.stats), 2)
    if (message) {
      message(paste0(
        "Marker call rate range in `geno.data` is: ",
        paste0(out$marker.callrate.stats, collapse = " ~ "),
        "."
      ))
    }
    out$maf.stats <- ASRgenomics:::maf(M = geno.data)
    out$maf.stats <- round(range(out$maf.stats), 2)
    if (message) {
      message(paste0(
        "Minor allele frequency in `geno.data` is: ",
        paste0(out$maf.stats, collapse = " ~ "),
        "."
      ))
    }
  }
  if (!is.null(Kinv)) {
    if (is.null(names(Kinv))) {
      names_of_Kinv <- paste0("Kinv", 1:length(Kinv))
    } else {
      names_of_Kinv <- names(Kinv)
    }
  }
  if (!is.null(Kinv)) {
    out$Kinv.condition <- lapply(Kinv, ASRgenomics:::Kinv.condition)
    out$Kinv.condition <- unlist(out$Kinv.condition)
    if (message) {
      message(paste0(
        "Condition of `Kinv` matrix(ces) is:",
        paste(
          "\n ",
          names_of_Kinv,
          ": ",
          out$Kinv.condition,
          sep = "",
          collapse = ", "
        ),
        "."
      ))
    }
  }
  return(out)
}
