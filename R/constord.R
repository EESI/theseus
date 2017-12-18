#' @import ggplot2
#' @importFrom phyloseq otu_table sample_data tax_table
NULL

#' Plots constrained ordination results
#'
#' Function 'conord' (constrained ordination) carries out constrained ordination
#' on a phyloseq object and plots the results as a ggplot2 object. Constrained
#' correspondence analysis ('CCA') and redundancy analysis ('RDA') are the two
#' methods currently implemented within this function.
#'
#' @param PS (required) A phyloseq object.
#' @param formula (required) Right-hand side of the model formula starting with
#'   a tilde ('~') and should not be placed in quotes. The current version
#'   requires at least two (2) constraining variables (see Details).
#' @param method Constrained ordination method to be applied. User may choose
#'   Constrained Correspondence Analysis (CCA) or Redundancy Analysis (RDA).
#'   Defaults to CCA.
#' @param facets Variable in sample_data(PS) to facet the plot by. Statement
#'   starts with a tilde ('~') and should not be placed in quotes.
#' @param scaling Scaling for species and site/sample scores in biplot. Options
#'   are the same as those found in the \code{\link[vegan]{scores}} function:
#'   "species" scaling (1) or 'site' scaling (2). The user should designate the
#'   appropriate scaling for thier intended analysis. Further information
#'   regarding scaling can be found in the Details below. Defaults to 2.
#' @param tax_level Taxonomic level to represent species composition using
#'   color. Defaults to 'Phylum'.
#' @param tax_n The number of taxonomic groups to identify using color (at the
#'   taxonomic level 'tax_level'). The most abundant tax_n will be selected. All
#'   other taxonomic groups will be collapsed into an additional 'Other'
#'   category for visualization. Defaults to 7.
#'
#' @details
#' \subsection{General}{
#' The current implementation of this function displays the first two
#'   constrained ordination axes. As such, the 'formula' argument must contain
#'   two or more constraining variables to return a valid plot; an error will be
#'   returned otherwise. Legendre and Legendre (1998, p. 587-592, 597-600,
#'   Table 11.1-11.5) provide thorough discussion on the constrained and
#'   unconstrained axes that result from constrained ordination methods.
#'   Selection of axes to plot, constrained or unconstrained, is a planned
#'   feature for a future release.
#' }
#'
#' \subsection{Comparison with phyloseq::plot_ordination}{
#' There are several differences between 'constord' and
#'   \code{\link[phyloseq]{plot_ordination}}, and they each have their own strengths. The highlights of 'constord' are:
#' \itemize{
#'   \item constraining variables are included in the 'constord' plot, a feature
#'     not currently present in phyloseq's approach, but still possible with
#'     some extra coding.
#'   \item 'constord' has argument 'scaling' which allows the user to select
#'     whether species scaling (1) or site scaling (2) is used when returning
#'     the scores to be plotting. Currently, *phyloseq::plot_ordination* returns
#'     site scaling (2). The choice of scaling is important and should be
#'     selected depending upon whether the goal is to compare the arrangement of
#'     sites or species (see Scalng section below).
#'   \item The aspect ratio of the ordination plots themselves are scaled
#'     according the ordination's eigenvalues to more accurately represent the
#'     distances between sites/samples, as described by Callahan et. al. (2016).
#'     Once again, this is easily included in
#'     \code{\link[phyloseq]{plot_ordination}}
#' }
#' }
#'
#' \subsection{Scaling}{
#' Species scaling (1) results in a distance biplot. The distance biplot is
#'   intended to enable the user to interpret the relationships between
#'   sites/samples. Site scaling (2) results in a correlation biplot. The
#'   correlation biplot enables the user to interpret the correlation between
#'   descriptors (species) within the ordination. Positions of sites/samples are
#'   not approximations of thier true locations; use species scaling (1) to
#'   interpret site/samples. A complete discussion of the implications of
#'   scaling (and interpretation of the ordination results) is provided in
#'   Legendre and Legendre (1998, p. 403-404, 585-587).
#' }
#'
#' @return A ggplot object.
#'
#' @references
#' Legendre, P. and Legendre, L. (1998) Numerical Ecology. 2nd English ed.
#'   Elsevier.
#' Callahan BJ, Sankaran K, Fukuyama JA et al. Bioconductor Workflow for
#'   Microbiome Data Analysis: from raw reads to community analyses
#'   [version 2; referees: 3 approved]. F1000Research 2016, 5:1492
#'   (doi: 10.12688/f1000research.8986.2)
#'
#' @seealso \code{\link[phyloseq]{ordinate}}
#'   \code{\link[phyloseq]{plot_ordination}} \code{\link[vegan]{cca}}
#'   \code{\link[vegan]{rda}} \code{\link[vegan]{scores}}
#'
#' @examples
#' \dontrun{
#' library(theseus)
#' data('WWTP_Impact')
#' p.co <- constord(PS=WWTP_Impact,
#'                  formula=~ log_NO3N + log_PO4,
#'                  method='RDA', facets=Position~Location, scaling=2)
#' p.co
#' }
#'
#' @export

constord <- function(PS,formula,method=c('CCA','RDA'),facets,scaling=2,tax_level='Phylum',tax_n=7){

  method <- match.arg(method)

  if (length(as.character(formula)) > 2)
    formula <- formula[-2]

  covariates <- all.vars(formula)

  if (!missing(facets)){
    facets_vars <- all.vars(facets)
    cov_subset <- c(covariates,facets_vars)
  }else{
    cov_subset <- covariates
  }


  sample_data(PS) <- sample_data(PS) %>%
    dplyr::select(cov_subset)

  sample_data(PS)$sample_id <- rownames(sample_data(PS))
  tax_table(PS) <- cbind(tax_table(PS),taxa_id=rownames(tax_table(PS)))

  PS <- phyloseq::prune_samples(sample_data(PS) %>%
                        tidyr::drop_na() %>%
                        rownames(.),PS)

  form <- stats::update.formula(formula, PS~.)

  ord <- phyloseq::ordinate(PS,method,formula=form)

  df <- phyloseq::psmelt(PS)
  top_taxa <- get_top_taxa(PS,tax_level,tax_n)
  df <- rename_column_to_other(df,tax_level,top_taxa)

  eig_p <- 100*ord$CCA$eig/sum(ord$CA$eig)
  aspect_ratio <- sqrt(ord$CCA$eig[2]/ord$CCA$eig[1])
  scores <- vegan::scores(ord,choices=seq_along(eig_p),scaling=scaling)

  sites <- data.frame(scores$sites)
  colnames(sites) <- paste0(colnames(sites),'_sites')
  sites$sample_id <- rownames(sites)
  df <- df %>% dplyr::left_join(sites,by='sample_id')

  species <- data.frame(scores$species)
  colnames(species) <- paste0(colnames(species),'_species')
  species$taxa_id <- rownames(species)
  df <- suppressWarnings(df %>% dplyr::left_join(species,by='taxa_id'))

  colnames(df) <- gsub(method,'axis',colnames(df))

  scaler <- df %>% dplyr::select(dplyr::starts_with('axis')) %>% unlist() %>% max(abs(.),na.rm=TRUE)
  bp <- data.frame(ord$CCA$biplot * scaler) %>%
    dplyr::mutate(covariate=rownames(.),
                  origin=0)
  colnames(bp) <- gsub(method,'axis',colnames(bp))

  p1 <- ggplot() +
    geom_vline(xintercept=0) +
    geom_hline(yintercept=0) +
    geom_point(data=df %>% dplyr::select_(~axis1_sites,~axis2_sites,~cov_subset) %>% dplyr::distinct(),
               aes_(x = ~axis1_sites, y = ~axis2_sites), alpha = 0.4,size=7) +
    geom_point(data=df %>% dplyr::select_(~axis1_species,~axis2_species,~Kingdom:Species) %>% dplyr::distinct(),
               aes_string(x = 'axis1_species', y = 'axis2_species', color = tax_level), alpha=0.7, size = 1) +
    geom_segment(data=bp,aes_(x=~origin,y=~origin,xend=~axis1,yend=~axis2),
                 arrow=grid::arrow(length=unit(0.03,'npc')),size=1.5,alpha=.5,color='darkred') +
    geom_text(data=bp,aes_(x=~axis1,y=~axis2,label=~covariate),fontface='bold') +
    guides(col = guide_legend(override.aes = list(size = 3))) +
    labs(x = sprintf("Axis1 [%s%% variance]", round(eig_p[1], 2)),
         y = sprintf("Axis2 [%s%% variance]", round(eig_p[2], 2))) +
    viridis::scale_color_viridis(discrete=TRUE) +
    theme_classic() +
    theme(panel.border = element_rect(color = "black", fill=NA, size=1),
          aspect.ratio = aspect_ratio) +
    ggtitle(method)

  if (!missing(facets))
    p1 <- p1 + facet_wrap(facets)

  p1

}