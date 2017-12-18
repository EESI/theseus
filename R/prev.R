#' @import ggplot2
#' @importFrom phyloseq otu_table sample_data tax_table
NULL

#' Create prevalence vs abundance plot
#'
#' Function 'prev' (prevalence plot) tabulates the prevalence and abundance of
#' each taxa in a phyloseq object and plots the results as a ggplot object.
#' This may assist the user in determining what filtering and preprocessing
#' steps should be taken regarding the removal of low count taxa.
#'
#' @param PS (required) A phyloseq object.
#' @param taxon Taxonomic level to be displayed. Defaults to 'Phylum'.
#' @param n_taxa The number of taxa to display at the given 'taxon' level.
#'   The most abundant taxa at that group are selected. Defaults to 10.
#' @param abund_threshold User selected value for the lowest acceptable total
#'   abundance. Defaults to 3.
#' @param prev_threshold User selected value for the lowest acceptable
#'   prevalence. Defaults to 5% of the number of samples in the phyloseq object,
#'   rounded up to the nearest integer.
#'
#' @details Low count taxa are often filtered from OTU tables to reduce possible
#'   error or noise. Examination of the raw (unfiltered) OTU table should be
#'   carried out to ensure that appropriate thresholds for prevalence (number of
#'   samples a taxa was observed in) and abundance (the total number of times a
#'   taxa was observed) are being selected. Function 'prev' plots each taxa
#'   according to thier prevalence and abundance within the dataset.
#'
#' @return A ggplot object.
#'
#' @references
#' Callahan BJ, Sankaran K, Fukuyama JA et al. Bioconductor Workflow for
#'   Microbiome Data Analysis: from raw reads to community analyses.
#'   F1000Research 2016, 5:1492 (doi: 10.12688/f1000research.8986.2)
#' Perraudeau F, Risso D, Street K et al. Bioconductor workflow for single-cell
#'   RNA sequencing: Normalization, dimensionality reduction, clustering, and
#'   lineage inference. F1000Research 2017, 6:1158
#'   (doi:10.12688/f1000research.12122.1)
#'
#' @examples
#' \dontrun{
#' library(theseus)
#' data('WWTP_Impact')
#' p.prev <- prev(WWTP_Impact, taxon="Phylum", n_taxa=10)
#' p.prev
#' }
#'
#' @export

prev <- function(PS,taxon,n_taxa=10,abund_threshold=3,prev_threshold=ceiling(0.05*phyloseq::nsamples(PS))){

  df <- data.frame(tax=phyloseq::taxa_names(PS),
                   abundance=phyloseq::taxa_sums(PS),
                   prevalence=colSums(phyloseq::otu_table(PS)>0)/phyloseq::nsamples(PS)) %>%
    dplyr::left_join(data.frame(tax=phyloseq::taxa_names(PS),phyloseq::tax_table(PS)),by='tax')

  if (!missing(taxon)){
    df <- df %>% dplyr::filter_(sprintf("%s %%in%% names(sort(table(.[,'%s']),decreasing=TRUE)[1:%s])",taxon,taxon,n_taxa))
  }

  p1 <- ggplot(df,aes_(~prevalence,~abundance)) +
    geom_point(alpha=.4) +
    scale_y_log10(breaks=10^(0:ceiling(log10(max(df$abundance))))) +
    theme(aspect.ratio=1,axis.text.x=element_text(angle=-90,vjust=.5,hjust=0)) +
    geom_hline(yintercept = abund_threshold, color='red', linetype=2) +
    geom_vline(xintercept = prev_threshold/phyloseq::nsamples(PS), color='red', linetype=2)

  if (!missing(taxon)){
    facet_form <- stats::as.formula(sprintf('~%s',taxon))
    p1 <- p1 + facet_wrap(facet_form)
  }

  p1

}

