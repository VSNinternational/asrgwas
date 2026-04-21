#' Generates phenotype-by-genotype plots for selected markers
#'
#' @description
#' Generates a plot, based on a given set of markers, for visualizing the
#' phenotype in relation to the underlying genotypes.
#' Both Gaussian and binomial data are supported.
#'
#' @param pheno.data A data frame with all relevant columns (factors and covariates)
#' and phenotypic responses to be used (default = \code{NULL}).
#' @param indiv A character with the name of the column in \code{pheno.data}
#' with the identification of treatments (genotypes) (default = \code{NULL}).
#' @param resp A character with the name of the numeric variable in \code{pheno.data} with the
#' response variable (default = \code{NULL}).
#' @param geno.data.sel A matrix with SNP data of form \eqn{n \times s},
#' with \eqn{n} individuals and \eqn{s} selected markers. Individual and marker names are
#' assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as: 0, 1, 2 (default = \code{NULL}).
#' @param type.plot A character indicating the type of plot. Options are: \code{"boxplot"} and
#' \code{"violin"} plot (default = \code{"boxplot"}).
#' @param family A character indicating the family of the distribution for \code{resp}.
#' Options are: \code{"gaussian"} and \code{"binomial"} (default = \code{"gaussian"}).
#' @param sample.label If \code{TRUE} the sample size of each genotypic state is added to the \eqn{x}
#' axis (default = \code{FALSE}).
#' @param sample.width If \code{TRUE} the size of the boxplots will be proportional to sample size
#' (default = \code{FALSE}).
#' @param maf.colour If \code{TRUE} the plot colour will be proportional to the minor allele
#' frequency (default = \code{FALSE}).
#' @param facet.dim A numeric vector indicating the number of rows and columns to be used
#' in the faceting process, \emph{e.g.}, \code{c(4, 1)} for 4 rows and 1 column (default = \code{NULL}).
#' @param scales A character indicating if the scales/levels in \eqn{x} and \eqn{y} axis should be the
#' same values across facets (panels). Options are: \code{"fixed"} \code{"free_x"}, \code{"free_y"}
#' and \code{"free"}. See \link[ggplot2]{facet_wrap} for more information (default = \code{NULL}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details
#'
#' These plots are useful to verify and explore the validity of reported associations.
#' For example: 1) a linear change is expected on significant markers for additive action,
#' 2) a non-linear pattern might indicate non-additive action,
#' 3) genotype classes might differ on their phenotypic distribution (\emph{e.g.}, heterogeneity of variances), and
#' 4) the influence of low MAF is clear and might indicate significant associations by chance.
#'
#' Check provided references for examples of these plots and their interpretation.
#'
#' @references
#' Croteau-Chonka D.C., Rogers A.J., Raj T., McGeachie M.J., Qiu W., et al. 2015.
#' Expression quantitative trait loci information improves predictive modeling of
#' disease relevance of non-coding genetic variation.
#' \emph{PLOS ONE} 10(10):e0140758.
#'
#' Galli G., Alves F.C., Morosini J.S., and Fritsche-Neto R. 2020.
#' On the usefulness of parental lines GWAS for predicting low heritability traits in tropical maize hybrids.
#' \emph{PLOS ONE} 15(2):e0228724.
#'
#' @return A ggplot object containing one or more phenotype-by-genotype plots.
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
#' # Get name and data of significant markers.
#' sign.markers <- gwasA$gwas.sel$marker
#' geno.data.sel <- gwas.data$geno.data[, sign.markers]
#'
#' # Simple marker plot call.
#' marker.plot(
#'  pheno.data = gwas.data$pheno.data, indiv = "Ind", resp = c("Sucrose"),
#'  geno.data.sel = geno.data.sel)
#'
#' # Manipulate some features (including a trait).
#' marker.plot(
#'  pheno.data = gwas.data$pheno.data, indiv = "Ind",
#'  resp = c("Sucrose", "Ethylene"), geno.data.sel = geno.data.sel,
#'  sample.label = TRUE, sample.width = TRUE, maf.colour = TRUE,
#'  scales = "free")
#' }
marker.plot <- function(
    pheno.data = NULL,
    indiv = NULL,
    resp = NULL,
    geno.data.sel = NULL,
    type.plot = c("boxplot", "violin"),
    family = c("gaussian", "binomial"),
    sample.label = FALSE,
    sample.width = FALSE,
    maf.colour = FALSE,
    facet.dim = NULL,
    scales = NULL,
    message = TRUE){
  .argument_class(.data = pheno.data, .class = "data.frame")
  pheno.data <- droplevels(pheno.data)
  .argument_class(.data = geno.data.sel, .class = c("matrix", "array"))
  if(is.null(rownames(geno.data.sel))){
    stop("Individual names not assigned to rows of `geno.data.sel`.")
  }
  if(is.null(colnames(geno.data.sel))){
    stop("Marker names not assigned to columns of `geno.data.sel`.")
  }
  .argument_class(.data = sample.label, .class = "logical")
  .argument_class(.data = sample.width, .class = "logical")
  .argument_class(.data = maf.colour, .class = "logical")
  .argument_class(.data = message, .class = "logical")
  .variable_class(
    .data = pheno.data, .mandatory =  TRUE, .variable = indiv, .class = "factor",
    .class.action = "message", .message = FALSE, .mutate = TRUE)
  .variable_class(
    .data = pheno.data, .mandatory =  TRUE, .variable = resp, .class = "numeric",
    .class.action = "message", .message = FALSE, .mutate = TRUE)
  type.plot <- match.arg(type.plot)
  family = match.arg(family)
  total.plots <- length(resp) * as.matrix(ncol(geno.data.sel))
  if (total.plots > 50){
    stop("Too many plots have been requested. Please select less than 50 marker-by-trait combinations to be plotted.")
  }
  geno.data.sel <- geno.data.sel[pheno.data[[indiv]], , drop = FALSE]
  melt.frame <- as.data.table(cbind.data.frame(pheno.data, geno.data.sel))
  melt.frame <- melt.data.table(
    melt.frame, id.vars = resp, measure.vars = colnames(geno.data.sel), na.rm = TRUE,
    variable.name = "Marker", value.name = "Genotype")
  if (sample.label) {
    merger.table <- aggregate(as.character(Genotype) ~ Marker + Genotype, data = melt.frame, length)
    names(merger.table)[ncol(merger.table)] <- "n"
    melt.frame <- merge.data.table(melt.frame, merger.table, by = c("Marker", "Genotype"))
    melt.frame$Genotype <- paste0(melt.frame$Genotype, '\n(', melt.frame$n, ")")
  }
  if (maf.colour) {
    merger.table <- data.frame(MAF = ASRgenomics:::maf(geno.data.sel))
    merger.table$Marker <- rownames(merger.table)
    melt.frame <- merge.data.table(melt.frame, merger.table, by = c("Marker"))
  }
  melt.frame <- melt.data.table(
    melt.frame, measure.vars = resp, na.rm = TRUE,
    variable.name = "Trait", value.name = "Phenotype")
  if (family == "binomial") {
    if (length(resp) > 1) {
      merger.table <- aggregate(as.character(Genotype) ~ Marker + Genotype + Phenotype, data = melt.frame, length)
      melt.frame <- merge.data.table(melt.frame, merger.table, by = c("Marker", "Genotype", "Phenotype"))
    }
    if (length(resp) == 1) {
      merger.table <- aggregate(as.character(Genotype) ~ Marker + Genotype + Trait + Phenotype, data = melt.frame, length)
      melt.frame <- merge.data.table(melt.frame, merger.table, by = c("Marker", "Genotype", "Trait",  "Phenotype"))
    }
    names(melt.frame)[ncol(melt.frame)] <- "n_by_phen"
  }
  melt.frame$Genotype <- as.factor(melt.frame$Genotype)
  facet.main <- "Marker"
  if (length(resp) > 1) {
    facet.secondary <- "Trait"
    scales.guess <- "free_y"
  } else {
    facet.secondary <- NULL
    scales.guess <- "fixed"
  }
  if (!is.null(facet.dim) & sample.label){
    if(facet.dim[1] > 1){
      scales.guess <- "free_x"
    }
    if(facet.dim[1] > 1 & length(resp) > 1){
      scales.guess <- "free"
    }
  }
  if (!is.null(facet.dim)){
    if (length(resp) > 1) {
      melt.frame$comp.facet <- paste0(melt.frame$Trait, " ", melt.frame$Marker)
      facet.main <- "comp.facet"
      facet.secondary <- NULL
      if (sample.label){
        scales.guess <- "free"
      } else {
        scales.guess <- "free_y"
      }
    } else {
      if (sample.label){
        scales.guess <- "free_x"
      }
    }
  }
  if (is.null(scales)){
    scales <- scales.guess
  }
  marker.plot <-
    ggplot(melt.frame, aes(x = Genotype, y = Phenotype)) +
    theme_light() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill = "#C0C0C0", colour = "#C0C0C0"))
  if (is.null(facet.dim)) {
    marker.plot <-
      marker.plot + facet_grid(reformulate(facet.main, facet.secondary), scales  = scales)
  } else {
    marker.plot <-
      marker.plot + facet_wrap(reformulate(facet.main, facet.secondary),
                               nrow = facet.dim[1], ncol = facet.dim[2], scales  = scales)
  }
  if (family == "binomial") {
    marker.plot <-
      marker.plot + geom_point(stat = "unique") +
      geom_vline(xintercept = levels(melt.frame$Genotype), alpha = 0.1, linetype = "dotted") +
      geom_text(label = melt.frame$n_by_phen, check_overlap = TRUE, size = 2.5,
                nudge_y = ifelse(melt.frame$Phenotype == 0, melt.frame$Phenotype - 0.1, melt.frame$Phenotype - 0.9)) +
      stat_summary(fun = "mean", geom = "point", size = 0.3) +
      scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75 , 1))
  }
  if (family == "gaussian") {
    if (type.plot == "boxplot")
      marker.plot <-
        marker.plot +
        geom_boxplot(outlier.fill = "white", outlier.size = .7,
                     outlier.shape = 21, alpha = 0.95,
                     varwidth = sample.width, fatten = 0.9, lwd = 0.3)
    if (type.plot == "violin")
      marker.plot <-
        marker.plot + geom_violin(alpha = 0.95) + stat_summary(fun = median, geom = "point", size = 0.3)
  }
  if (maf.colour){
    marker.plot <-
      marker.plot + aes(colour = MAF) +
      scale_colour_viridis_c(direction = -1)
  }
  return(marker.plot)
}
