#' Draws the genetic map
#'
#' @description Generates the genetic map for visualization.
#'
#' @param map.data A data frame containing the map information (\code{"marker"}, \code{"chrom"},
#'  and \code{"pos"}) (default = \code{NULL}).
#' @param tag.markers A (non-mandatory) vector with the identification (name) of markers to be highlighted
#' in the map (default = \code{NULL}).
#' @param tag.repel If \code{TRUE} tags are automatically repelled (default = \code{TRUE}).
#' @param tag.colour A character with a colour name (\emph{e.g.}, \code{"black"})
#' to be used in the colouring of selected markers labels (default = \code{"black"}).
#' @param tag.size A numeric value (>= 0) indicating label sizes for
#'  selected marker names (default = \code{1.3}).
#' @param marker.colour A character with a colour name (\emph{e.g.}, \code{"blue"})
#' to be used in marker colouring (default = \code{"#105E26"}).
#' @param marker.width A numeric value (>= 0) indicating the marker width (represented by a line;
#' default = \code{0.05}).
#' @param marker.alpha A numeric value indicating the markers' transparency (default = \code{0.08}).
#' @param chrom.colour A character with a colour name (\emph{e.g.}, \code{"white"})
#' to be used in inner chromosome colouring (default = \code{"#DEDEDE"}).
#' @param chrom.width A numeric value (> 0) indicating the chromosome width (default = \code{4}).
#' @param chrom.contour A character with a colour name (\emph{e.g.} \code{"#DEDEDE"})
#' to be used in outer chromosome colouring (default = \code{"black"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A ggplot object containing the genetic map.
#'
#' @details This function was developed in collaboration with Khaled Al-Shamaa (ICARDA/CGIAR).
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
#' # Get genetic map with the significant markers.
#' map.plot(map.data = gwasA$gwas.all,
#'   tag.markers = gwasA$gwas.sel$marker)
#'
#' # Customizing a bit.
#' map.plot(map.data = gwasA$gwas.all,
#'   tag.markers = gwasA$gwas.sel$marker,
#'   marker.colour = "white", marker.width = 0.2,
#'   marker.alpha = 0.1, chrom.colour = "grey",
#'   chrom.contour = "grey", tag.colour = "#0072B2")
#' }
map.plot <- function(
    map.data = NULL,
    tag.markers = NULL,
    tag.repel = TRUE,
    tag.colour = "black",
    tag.size = 1.8,
    marker.colour = "#105E26",
    marker.width = 0.05,
    marker.alpha = 0.08,
    chrom.colour = "#DEDEDE",
    chrom.width = 4,
    chrom.contour = "#DEDEDE",
    message = TRUE
){
  .argument_class(.data = map.data, .class = "data.frame")
  if (!is.null(tag.markers)){
    .argument_class(.data = tag.markers, .class = "character")
  }
  .argument_class(.data = marker.colour, .class = "character")
  .argument_class(.data = chrom.colour, .class = "character")
  .argument_class(.data = chrom.contour, .class = "character")
  .argument_class(.data = tag.colour, .class = "character")
  map.data <- droplevels(map.data)
  mandatory.variables <- c("marker", "chrom", "pos")
  if (!all(mandatory.variables %in% colnames(map.data))){
    stop(paste0("The provided `map.data` is missing one or of the following variables: ",
                paste0(mandatory.variables, collapse = ", ")))
  }
  if (marker.width < 0 | !is.numeric(marker.width)){
    stop("The `marker.width` argument should be a numeric value >= 0.")
  }
  if (marker.alpha < 0 | marker.alpha > 1 | !is.numeric(marker.alpha)){
    stop("The `marker.alpha` argument should be a numeric value following the condition 0 <= marker.alpha <= 1.")
  }
  if (chrom.width <= 0 | !is.numeric(chrom.width)){
    stop("The `chrom.width` argument should be a numeric value > 0.")
  }
  if (tag.size < 0 | !is.numeric(tag.size)){
    stop("The `tag.size` argument should be >= 0.")
  }
  .argument_class(.data = tag.repel, .class = "logical")
  .argument_class(.data = message, .class = "logical")
  map.data$chrom <- as.factor(map.data$chrom)
  map.data$chrom_pos <- as.numeric(factor(map.data$chrom))
  if (!is.null(tag.markers)){
    levels(map.data$chrom) <- c(levels(map.data$chrom), "")
  }
  map.data.max = aggregate(pos ~ chrom, map.data, max)
  map.data.sel = subset(map.data, marker %in% tag.markers)
  base.plot <- ggplot(map.data, aes(chrom, pos))
  base.plot <- base.plot +
    geom_segment(
      data = map.data.max,
      mapping =
        aes(x = chrom,
            xend = chrom,
            y = 0,
            yend = pos),
      linewidth = chrom.width,
      color = chrom.contour,
      lineend = "round"
    )
  base.plot <- base.plot +
    geom_segment(
      data = map.data.max,
      mapping =
        aes(x = chrom,
            xend = chrom,
            y = 0,
            yend = pos),
      linewidth = chrom.width - 0.02,
      color = chrom.colour,
      lineend = "round"
    )
  base.plot <- base.plot +
    geom_segment(
      data = map.data,
      mapping =
        aes(x = chrom_pos - marker.width,
            xend = chrom_pos + marker.width,
            y = pos,
            yend = pos),
      colour = marker.colour,
      alpha = marker.alpha
    )
  base.plot <- base.plot +
    geom_point(
      data = map.data.sel,
      color = tag.colour,
      size = 6, shape = "_"
    )
  if (tag.repel){
    base.plot <- base.plot +
      geom_text_repel(
        data = map.data.sel,
        aes(chrom, pos, label = marker),
        hjust = 0,
        color = tag.colour,
        nudge_x = 0.2, direction = "y",
        size = tag.size
      )
  }
  if (!tag.repel){
    base.plot <- base.plot +
      geom_text(
        data = map.data.sel,
        aes(chrom, pos, label = marker),
        hjust = 0,
        color = tag.colour,
        nudge_x = 0.2,
        size = tag.size
        )
  }
  base.plot <- base.plot +
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.ticks.x = element_blank()) +
    scale_x_discrete(
      name = "Chromosome / Linkage group",
      drop = FALSE,
      ) +
    scale_y_reverse(
      name = "Position",
      labels = function(x) format(x, scientific = TRUE))
  return(base.plot)
}
