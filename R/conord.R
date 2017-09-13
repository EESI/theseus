#' @import ggplot2
NULL

#' Plots quality assessment
#'
#' Description paragraph
#'
#' @param PS (required) asdf
#' @param formula (required) asdf
#' @param method asdf. Defaults to CCA.
#' @param facets asdf
#' @param scaling asdf. Defaults to 2.
#' @param tax_level asdf. Defaults to Phylum.
#' @param tax_n asdf. Defaults to 7.
#'
#' @return A ggplot object.
#'
#' @references
#' X
#' Y
#' Z
#'
#' @seealso \code{\link[phyloseq]{ordinate}} \code{\link[vegan]{cca}} \code{\link[vegan]{rda}}
#'
#' @examples
#' \dontrun{
#' conord(PS,~ SL_NPOC + SL_NO3 + SL_NH4,facets=~reactor.cat,tax_level='Class',scaling=2,method='RDA')
#' }
#'
#' @export

conord <- function(PS,formula,method=c('CCA','RDA'),facets,scaling=2,tax_level='Phylum',tax_n=7){

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

  PS <- prune_samples(sample_data(PS) %>%
                        tidyr::drop_na() %>%
                        rownames(.),PS)

  form <- update.formula(formula, PS~.)

  ord <- ordinate(PS,method,formula=form)

  df <- psmelt(PS)
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
    geom_point(data=df %>% dplyr::select(axis1_sites,axis2_sites,cov_subset) %>% dplyr::distinct(),
               aes(x = axis1_sites, y = axis2_sites), alpha = 0.4,size=7) +
    geom_point(data=df %>% dplyr::select(axis1_species,axis2_species,Kingdom:Species) %>% dplyr::distinct(),
               aes_string(x = 'axis1_species', y = 'axis2_species', color = tax_level), alpha=0.7, size = 1) +
    geom_segment(data=bp,aes(x=origin,y=origin,xend=axis1,yend=axis2),
                 arrow=grid::arrow(length=unit(0.03,'npc')),size=1.5,alpha=.5,color='darkred') +
    geom_text(data=bp,aes(x=axis1,y=axis2,label=covariate),fontface='bold') +
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