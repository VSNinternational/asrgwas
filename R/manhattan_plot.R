#' Draws one or more Manhattan plots for GWAS marker data results
#'
#' @description Generates one or more Manhattan plots for visualization.
#'
#' @param gwas.table A data frame of one or more (vertically stacked) GWAS analysis with (as a minimum) the
#' following variables: \code{"marker"}, \code{"chrom"}, \code{"pos"},
#' and \code{"p.value"} of markers (default = \code{NULL}).
#' @param pvalue.thr A numeric value with the \code{p-value} threshold to identify significant markers (\emph{e.g.}, \code{5e-4}).
#' The value will be transformed using \code{-log10(pvalue.thr)} for plotting (default = \code{NULL}).
#' @param tag.table A data frame of tags (annotations) for scores. The columns \code{"marker"} and
#' \code{"tag"} are required. The column \code{"tag"} should contain the
#' text to be annotated in each marker score (\emph{e.g.}, "Candidate region 3"). Columns
#' \code{"nudge.x"} and \code{"nudge.y"} are optional to inform position nudges
#' on the text and should be on the base pair number and \code{log(p-value)} scale,
#' respectively (default = \code{NULL}).
#' @param tag.repel If \code{TRUE}, \code{nudge.x} and \code{nudge.y} from \code{tag.table} are ignored
#' and tags are automatically repelled (default = \code{TRUE}).
#' @param point.colour A character with a colour name (\emph{e.g.}, \code{"blue"}) or a column name
#' from \code{gwas.table} to be used as the index for point colouring (default = \code{"chrom"}).
#' @param point.size A numeric value (greater than 0) indicating the point size
#' (default = \code{1}).
#' @param point.alpha A numeric value indicating the points' transparency (default = \code{1}).
#' @param padding A numeric value that will be used as the distance between chromosomes (default = \code{0}).
#' @param collate If \code{TRUE} the distance between the last marker of a chromosome and
#' the first marker of the next chromosome is set to 1 base pair (default = \code{TRUE}).
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
#' @param scales A character indicating if the scales/levels in \eqn{x} and \eqn{y} axis should be the
#' same values across facets (panels). Options are: \code{"fixed"} \code{"free_x"}, \code{"free_y"},
#' and \code{"free"}. See \link[ggplot2]{facet_wrap} for more information
#' (default = \code{"free_y"}).
#' @param legend.position A character indicating the position for the legend in the displayed plot.
#' Options are: \code{"bottom"}, \code{"top"}, \code{"left"}, \code{"right"}, \code{"none"}, or the
#' corresponding vector of coordinates (such as: \code{c(0.95, 0.05)} for \eqn{x} and \eqn{y}, respectively)
#' (default = \code{"none"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A ggplot object containing one or more Manhattan plots.
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
#'   pvalue.thr = 0.0005, bonferroni = FALSE, workspace = "2Gb")
#'
#' # Simple Manhattan plot call.
#' manhattan.plot(gwas.table = gwasA$gwas.all)
#'
#' # Manipulate some features.
#' manhattan.plot(
#'   gwas.table = gwasA$gwas.all,
#'   pvalue.thr = 0.0005, point.alpha = 0.80, padding = 3e6)
#'
#' # Adding tags.
#' gwasA$gwas.sel$tag <-
#'   c("Putative gene 1", "Putative gene 2",
#'     "Putative gene 3", "Putative gene 4", "Putative gene 5")
#'
#' manhattan.plot(
#'   gwas.table = gwasA$gwas.all,
#'   pvalue.thr = 0.0005, point.alpha = 0.80, padding = 3e6,
#'   tag.table = gwasA$gwas.sel, tag.repel = TRUE)
#' }
manhattan.plot <- function(
    gwas.table = NULL,
    pvalue.thr = NULL,
    tag.table = NULL,
    tag.repel = TRUE,
    point.colour = "chrom",
    point.size = 1,
    point.alpha = 1,
    collate = TRUE,
    padding = 0,
    gwas.index = NULL,
    facet.main = NULL,
    facet.secondary = NULL,
    facet.dim = NULL,
    scales = "free_y",
    legend.position = "none",
    message = TRUE
){
  if (!is.character(point.colour)){
    stop("The `point.colour` argument should be of class character.")
  }
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
  mandatory.variables <- c("marker", "chrom", "pos", "p.value",
                           gwas.index, facet.main, facet.secondary)
  colour.aes <- point.colour %in% colnames(gwas.table)
  if (colour.aes){
    mandatory.variables <- append(mandatory.variables, point.colour)
  }
  if (!all(mandatory.variables %in% colnames(gwas.table))){
    stop(paste0("The provided `gwas.table` is missing one or of the following variables: ",
                paste0(mandatory.variables, collapse = ", ")))
  }
  if (!is.null(tag.table)){
    .argument_class(.data = tag.table, .class = "data.frame")
    mandatory.tag.variables <- c("marker", "tag", facet.main, facet.secondary)
    if (!all(mandatory.tag.variables %in% colnames(tag.table))){
      stop(paste0("The provided `tag.table` is missing one or more mandatory variables: ",
                  paste0(mandatory.tag.variables, collapse = ", ")))
    }
  }
  if (!is.null(pvalue.thr)){
    if (!is.numeric(pvalue.thr)){
      stop("The value(s) provided in the `pvalue.thr` argument should be of class `numeric`.")
    }
  }
  if (point.size < 0){
    stop("The `point.size` argument should follow the condition 0 <= point.size <= 1.")
  }
  if (point.alpha <= 0 || point.alpha > 1){
    stop("The `point.alpha` argument should follow the condition 0 < point.size <= 1.")
  }
  if (!is.null(facet.dim)){
    if (length(facet.dim) != 2 || !is.numeric(facet.dim)){
      stop("The `facet.dim` aregument should be a vector with two numeric values (number of rows and columns).")
    }
  }
  .argument_class(.data = tag.repel, .class = "logical")
  .argument_class(.data = message, .class = "logical")
  if (!is.data.table(gwas.table)) {
    gwas.table <- as.data.table(gwas.table)
  }
  if (!is.null(tag.table)){
    if (!is.data.table(tag.table)) tag.table <- as.data.table(tag.table)
    if (!"nudge.x" %in% colnames(tag.table)) {tag.table[, nudge.x := 0]}
    if (!"nudge.y" %in% colnames(tag.table)) {tag.table[, nudge.y := 0]}
  }
  if (!is.null(gwas.index)){
    n.frames <- length(unique(gwas.table[[gwas.index]]))
  } else {
    n.frames = 1
    gwas.table$gwas.index <- "gwas001"
  }
  if (message){
    message(paste0("A total of ", n.frames, " analysis/analyses reported in `gwas.table`."))
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
    if (facet.dim[1] *  facet.dim[2] != n.panels & message){
      message("The number of rows and columns provided in `facet.dim` does not correspond to the number of facets provided.")
    }
  }
  if (n.frames != n.panels & !colour.aes  & message){
    message("The use of a single colour is not recommended when more than one analysis is presented by panel.")
  }
  if (n.frames < length(pvalue.thr)){
    stop("The number of thresholds is larger than the number of GWAS analyses provided.")
  }
  if ((!is.null(pvalue.thr)) & n.frames > length(pvalue.thr) & message){
    message(paste0("The number of thresholds is lower than number of analysis, they will be recycled in the provided order."))
  }
  gwas.table[, x:= marker.position(
    map.data = .SD[, c("pos", "chrom")], chrom = "chrom", pos = "pos",
    collate = collate, padding = padding), by = gwas.index]
  chrom.breaks <- tapply(X = gwas.table[["x"]], INDEX = gwas.table[["chrom"]],
                         FUN = function(x) (max(x) - min(x))/2 + min(x))
  tmp.row.count <- nrow(gwas.table)
  if (!is.null(tag.table)){
    gwas.table <- merge.data.table(tag.table, gwas.table, all = TRUE)
    if(tmp.row.count != nrow(gwas.table)){
      stop("Something went wrong merging `gwas.table` and `tag.table`. It is possible that the common columns are different.")
    }
  }
  if (!is.null(pvalue.thr)){
    pvalue.thr <- -log10(pvalue.thr)
    if (any(!is.null(facet.main), !is.null(facet.secondary))){
      unique.groups <- gwas.table[,.GRP, by = c(facet.main, facet.secondary)]
      pvalue.thr <- cbind(unique.groups, "pvalue.thr" = pvalue.thr)
      pvalue.thr[, "GRP" := NULL]
      gwas.table <- merge.data.table(gwas.table, pvalue.thr, all.x = TRUE, by = c(facet.main, facet.secondary))
    } else{
      gwas.table <- cbind(gwas.table, "pvalue.thr" = pvalue.thr)
    }
    if(tmp.row.count != nrow(gwas.table)){
      stop("Something went wrong merging `gwas.table` and `pvalue.thr`. It is possible that the common columns are different.")
    }
  }
  if (colour.aes){
  base.plot <-
    ggplot(data = gwas.table,
           mapping = aes(x = x, y = -log10(p.value),
                         colour = factor(!!as.name(point.colour)))) +
    labs(colour = point.colour)
  } else {
    base.plot <-
      ggplot(data = gwas.table,
             mapping = aes(x = x, y = -log10(p.value)))
  }
  base.plot <- base.plot +
    theme_light() +
    scale_colour_viridis_d() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = legend.position) +
    scale_x_continuous(name = "Chromosome / Linkage group",
                       breaks = chrom.breaks, labels = names(chrom.breaks)) +
    scale_y_continuous(name = "-log10(p-value)")
  if(!is.null(pvalue.thr)){
  base.plot <-
    base.plot +
    geom_line(aes(y = gwas.table[["pvalue.thr"]]), colour = "black", linetype = "dashed", alpha = 0.7, linewidth = 0.4)
  }
  if(!is.null(tag.table)){
    if (!tag.repel){
      base.plot <-
        base.plot +
        geom_text(data = na.omit(gwas.table), mapping =  aes(x = x, y = -log10(p.value)),
                  label = na.omit(gwas.table[["tag"]]), colour = "black",
                  nudge_x = na.omit(gwas.table[["nudge.x"]]), nudge_y = na.omit(gwas.table[["nudge.y"]]))
    }
    if (tag.repel) {
      base.plot <-
        base.plot +
        geom_text_repel(data = na.omit(gwas.table), mapping =  aes(x = x, y = -log10(p.value)),
                  label = na.omit(gwas.table[["tag"]]), colour = "black")
      }
    }
  if (colour.aes){
    base.plot <-
      base.plot +
      geom_point(size = point.size, alpha = point.alpha, shape = 15)
  } else {
    base.plot <-
      base.plot +
      geom_point(size = point.size, alpha = point.alpha, colour = point.colour, shape = 15)
  }
  if (any(!is.null(facet.main), !is.null(facet.secondary))){
    if (is.null(facet.dim)) {
      base.plot <-
        base.plot + facet_grid(reformulate(facet.main, facet.secondary), scales = scales)
    } else {
      base.plot <-
        base.plot + facet_wrap(reformulate(facet.main, facet.secondary),
                               nrow = facet.dim[1], ncol = facet.dim[2], scales = scales)
    }
  }
  return(base.plot)
}
