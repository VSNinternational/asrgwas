#' Generates continuous relative position of markers
#'
#' @param map.data A data frame with each marker's name, chromosome, and position.
#' Variable names \strong{must} be "marker", "chrom", and "pos" (default = \code{NULL}).
#' @param collate If \code{TRUE} the chromosomes will be separated by one base pair (default = \code{TRUE}).
#' @param padding A \code{numeric} value that will be used as the distance between chromosomes (default = \code{0}).
#'
#' @keywords internal
#'
#' @return A vector with the continuous relative positions of markers.
#'
#' @examples
#' map.data <- ASRgenomics:::dummy.map_(marker.id = 1:20)
#' map.data$chrom[11:20] <- 2
#' map.data$pos <- ASRgwas:::marker.position(
#'  map.data, collate = TRUE, padding = 5)
#' map.data
marker.position <- function(map.data = NULL,
                            chrom = "chrom", pos = "pos",
                            collate = TRUE, padding = 0,
                            message = TRUE){
  .argument_class(.data = map.data, .class = "data.frame")
  if (length(unique(map.data[[chrom]])) == 1){
    if (message) {
      message("The map only contains one chromosome, no modifications done to positions.")
    }
    return(map.data[[pos]])
  }
  .argument_class(.data = collate, .class = "logical")
  if(!is.numeric(padding) || padding < 0) {
    stop("The `padding` argument should be a numeric value larger than 0.")
  }
  if (collate){
    chrom.min <- tapply(X = map.data[[pos]], INDEX = map.data[[chrom]], FUN = min) - 1
    chrom.min[1] <- 0
    chrom.min <- cumsum(chrom.min)
  } else {
    chrom.min <- tapply(X = rep(0, length(map.data[[chrom]])), INDEX = map.data[[chrom]], sum)
  }
  chrom.max <- tapply(X = map.data[[pos]], INDEX = map.data[[chrom]], FUN = max) + padding
  chrom.max.names <- names(chrom.max)
  chrom.max <- shift(cumsum(chrom.max), fill = 0)
  names(chrom.max) <- chrom.max.names
  pos <- map.data[[pos]] + chrom.max[map.data[[chrom]]] - chrom.min[map.data[[chrom]]]
  return(pos)
}
