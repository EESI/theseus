#' @import ggplot2
#' @importFrom phyloseq otu_table sample_data tax_table
NULL

#' Plots constrained ordination results
#'
#' Function conord performs Constrained Correspondence Analysis (CCA) or Redundancy Analysis (RDA) and plots the results as a ggplot object.
#'
#' @param PS (required) A phyloseq object.
#' @param taxon asdf
#' @param n_taxa adsf
#' @param theshold asdf
#'
#' @return A ggplot object.
#'
#' @references
#' Perraudeau F, Risso D, Street K et al. Bioconductor workflow for single-cell RNA sequencing:
#' Normalization, dimensionality reduction, clustering, and lineage inference.
#' F1000Research 2017, 6:1158 (doi: 10.12688/f1000research.12122.1)
#'
#'
#' @examples
#' \dontrun{
#' sometest
#' }
#'
#' @export

prev <- function(PS,taxon='Phylum',n_taxa=10,threshold=3){

  facet_form <- stats::as.formula(sprintf('~%s',taxon))

  df <- data.frame(tax=phyloseq::taxa_names(PS),
                   abundance=phyloseq::taxa_sums(PS),
                   prevalence=colSums(otu_table(PS)>0)/phyloseq::nsamples(PS)) %>%
    dplyr::left_join(data.frame(tax=phyloseq::taxa_names(PS),tax_table(PS)),by='tax') %>%
    dplyr::filter_(sprintf("%s %%in%% names(sort(table(.[,'%s']),decreasing=TRUE)[1:%s])",taxon,taxon,n_taxa))

  ggplot(df,aes_(~abundance,~prevalence)) +
    geom_hline(yintercept=threshold/phyloseq::nsamples(PS),color='red',linetype=2) +
    geom_point(alpha=.4) +
    facet_wrap(facet_form) +
    scale_x_log10(breaks=10^(0:ceiling(log10(max(df$abundance))))) +
    theme(aspect.ratio=.5,axis.text.x=element_text(angle=-90,vjust=.5,hjust=0))

}

