#' Draws one or more Q-Q plots for GWAS marker data results
#'
#' @description Generates one or more quantile-quantile (Q-Q) plots for visualization.
#'
#' @param gwas.table A data frame of one or more (vertically stacked) GWAS analyses with (as a minimum)
#' the \code{"p.value"} column of marker results (default = \code{NULL}).
#' @param point.colour A character with a colour name (\emph{e.g.}, \code{"blue"}) or a column name
#' from \code{gwas.table} to be used as the index for point colouring (default = \code{"#0072B2"}).
#' @param point.size A numeric value (greater than 0) indicating the point size
#' (default = \code{1}).
#' @param point.alpha A numeric value indicating the points' transparency (default = \code{1}).
#' @param trend.colour A character with a colour name (\emph{e.g.}, \code{"blue"}) to be used
#' for colouring of the central line (expected value) (default = \code{"black"}).
#' @param coord.equal If \code{TRUE} the \eqn{x} and \eqn{y} coordinates will have the same limits
#' and a ratio of 1 (default = \code{TRUE}).
#' @param gwas.index A character indicating a column name in \code{gwas.table} identifying each GWAS analysis.
#' This is only necessary when more than one GWAS analysis is present in \code{gwas.table} (default = \code{NULL}).
#' @param facet.main A character indicating a column name in \code{gwas.table} to be used
#' for the main faceting of panels (if more than one GWAS analysis is present) (default = \code{NULL}).
#' @param facet.secondary A character indicating a column name in \code{gwas.table} to be used
#' as secondary faceting of panels (if more than one GWAS analysis is present) (default = \code{NULL}).
#' @param facet.dim A numeric vector indicating the number of rows and columns to be used
#' in the faceting process, \emph{e.g.}, \code{c(4, 1)}, for 4 rows and 1 column. This argument only works
#' if \code{facet.main} is provided, and works best if \code{facet.secondary} is \code{NULL}
#' (default = \code{NULL}).
#' @param legend.position A character indicating the position for the legend in the displayed plot.
#' Options are: \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"none"}, or the
#' corresponding vector of coordinates (such as \code{c(0.95, 0.05)} for \eqn{x} and \eqn{y}, respectively)
#' (default = \code{"none"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A ggplot objects containing one or more Q-Q plots.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Prepare Apricot dataset.
#' gwas.data <- pre.gwas(
#'   pheno.data = pheno.apricot, indiv = "Ind", resp = "Sucrose",
#'   geno.data = geno.apricot, map.data = map.apricot, maf = 0.05,
#'   marker.callrate = 0.40, impute = TRUE)
#'
#' # Run GWAS.
#' gwasA <- gwas.asreml(
#'   pheno.data = gwas.data$pheno.data, resp = "Sucrose", gen = "Ind",
#'   fixedf = "Lots", residual = "Lots",
#'   Kinv = gwas.data$Kinv, Q = gwas.data$Q, npc = 3,
#'   geno.data = gwas.data$geno.data, map.data = gwas.data$map.data,
#'   pvalue.thr = 0.0005, workspace = "2Gb")
#'
#' # Simple QQ plot call.
#' qq.plot(gwas.table = gwasA$gwas.all)
#'
#' # Manipulate some features.
#' qq.plot(
#'   gwas.table = gwasA$gwas.all, point.colour = "green",
#'   point.alpha = 0.5, point.size = 2)
#' }
qq.plot <- function(
    gwas.table = NULL,
    point.colour = "#0072B2",
    point.size = 1,
    point.alpha = 1,
    trend.colour = "grey",
    coord.equal = TRUE,
    gwas.index = NULL,
    facet.main = NULL,
    facet.secondary = NULL,
    facet.dim = NULL,
    legend.position = "none",
    message = TRUE
){
  .argument_class(.data = point.colour, .class = "character")
  .argument_class(.data = trend.colour, .class = "character")
  if (!is.null(facet.main)){
    if (!is.character(facet.main)){
      stop("The `facet.main` argument should be of class character.")
    }
  }
  if (!is.null(facet.secondary)){
    if (!is.character(facet.secondary)){
      stop("The `facet.secondary` argument should be of class character.")
    }
  }
  .argument_class(.data = gwas.table, .class = "data.frame")
  gwas.table <- droplevels(gwas.table)
  mandatory.variables <- c("p.value", facet.main, facet.secondary)
  colour.aes <- point.colour %in% colnames(gwas.table)
  if (colour.aes){
    mandatory.variables <- append(mandatory.variables, point.colour)
  }
  if (!all(mandatory.variables %in% colnames(gwas.table))){
    stop(paste0("The provided `gwas.table` is missing one or of the following variables: ",
                paste0(mandatory.variables, collapse = ", ")))
  }
  if (point.size <= 0){
    stop("The `point.size` argument should follow the condition 0 < point.size <= 1.")
  }
  if (point.alpha <= 0 || point.alpha > 1){
    stop("The `point.alpha` argument should follow the condition 0 < point.size <= 1.")
  }
  if (!is.null(facet.dim)){
    if (length(facet.dim) != 2 || !is.numeric(facet.dim)){
      stop("The `facet.dim` argument should be a vector with two numeric values (number of rows and columns).")
    }
  }
  .argument_class(.data = coord.equal, .class = "logical")
  .argument_class(.data = message, .class = "logical")
  if (!is.data.table(gwas.table)) {
    gwas.table <- as.data.table(gwas.table)
  }
  if (!is.null(gwas.index)){
    n.frames <- length(unique(gwas.table[[gwas.index]]))
  } else {
    n.frames = 1
    gwas.table$gwas.index <- "gwas001"
  }
  if (message){
    message(paste0("A total of ", n.frames, " analysis/analyses included in `gwas.table`."))
  }
  if (n.frames == 1 & max(table(gwas.table[["marker"]])) > 1){
    if (message){
      message("It is possible that more than one analysis is reported in `gwas.table` but no `gwas.index` argument was provided.")
    }
  }
  if (any(!is.null(facet.main), !is.null(facet.secondary))){
    n.panels <- nrow(gwas.table[,.GRP, by = c(facet.main, facet.secondary)])
  } else {
    n.panels <- 1
  }
  if (!is.null(facet.dim)){
    if (facet.dim[1] *  facet.dim[2] != n.panels){
      warning("The number of rows and columns provided in `facet.dim` does not correspond to the number of facets provided.")
    }
  }
  if (!colour.aes & n.frames != n.panels & message){
    message("The use of a single colour is not recommended when more than one analysis is presented by panel.")
  }
  gwas.table <- gwas.table[order(p.value), .SD, by = gwas.index]
  gwas.table[, obs := -log10(p.value)]
  gwas.table[, exp := -log10(1:nrow(.SD)/(nrow(.SD) + 1)), by = gwas.index]
  if (colour.aes){
    base.plot <-
      ggplot(data = gwas.table,
             mapping = aes(x = exp, y = obs, colour = factor(!!as.name(point.colour)))) +
      labs(colour = point.colour, x = 'Theoretical  -log10(p-value)', y = 'Observed  -log10(p-value)')
  } else {
    base.plot <-
      ggplot(data = gwas.table, mapping = aes(x = exp, y = obs)) +
      labs(x = 'Theoretical  -log10(p-value)', y = 'Observed  -log10(p-value)')
  }
  base.plot <- base.plot +
    theme_light() +
    scale_colour_viridis_d() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = legend.position)
  base.plot <- base.plot +
    geom_abline(
      slope = 1, intercept = 0, linetype = "dashed", linewidth = .8, colour = trend.colour)
  if (colour.aes){
    base.plot <- base.plot +
      geom_point(size = point.size, alpha = point.alpha, shape = 16)
  } else {
    base.plot <- base.plot +
      geom_point(size = point.size, alpha = point.alpha, shape = 16,
                 colour = point.colour)
  }
  if (coord.equal) {
    xy.lim <- max(c(layer_scales(base.plot)$x$range$range[2] ,
                    layer_scales(base.plot)$y$range$range[2]))
    base.plot <- base.plot + coord_equal(xlim = c(0, xy.lim), ylim = c(0, xy.lim))
  }
  if (any(!is.null(facet.main), !is.null(facet.secondary))){
    if (is.null(facet.dim)) {
      base.plot <-
        base.plot + facet_grid(reformulate(facet.main, facet.secondary))
    } else {
      base.plot <-
        base.plot + facet_wrap(reformulate(facet.main, facet.secondary),
                               nrow = facet.dim[1], ncol = facet.dim[2])
    }
  }
  return(base.plot)
}
