#' Freshwater stream microbiome data
#'
#' WWTP_Impact is a phyloseq object containing 16S rDNA amplicon sequencing data
#' from samples collected at 6 sites on 2 days (12 total samples) along
#' Wissahickon Creek and Sandy Run in southeastern Pennsylvania, USA. Sequencing
#' data was processed as described in Price et. al. (2017), and combined with
#' chemical data to create a phyloseq object. The raw data and scripts used to
#' generate this phyloseq object (as well as the phyloseq object itself) can be
#' obtained from the author's GitHub repository (see url below).
#'
#' @docType data
#' @name WWTP_Impact
#' @usage WWTP_Impact
#' @keywords datasets
#'
#' @format A phyloseq object containing sample_data, otu_table, tax_table, and
#'   phy_tree slots.
#'
#' @references Price, J. R., Ledford, S. H., Ryan, M. O., Toran, L., Sales,
#'    C. M. (2017). Wastewater treatment plant effluent introduces recoverable
#'    shifts in microbial community composition in receiving streams. Sci. Total
#'    Environ. doi:10.1016/j.scitotenv.2017.09.162.
#'
#' @source \url{https://github.com/JacobRPrice/WWTP_Impact_on_Stream}
NULL