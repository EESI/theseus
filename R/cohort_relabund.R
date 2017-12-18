#' @import ggplot2
#' @importFrom phyloseq transform_sample_counts plot_bar
NULL

#' Creates a relative abundance cohort plot
#'
#' This function plots the relative abundance of taxa within a phyloseq
#' object 'PS' according to thier pre-determined cohort memberships.
#' See Details for more information.
#'
#' @param PS (required) A phyloseq object.
#' @param xvar (required) Variable in sample_data(PS) to be displayed on the
#'   x-axis. Defaults to 'SampleID'.
#' @param taxfill Taxonomic level to display in plot. Defaults to 'Phylum'.
#' @param comp1 (required) First comparison (hence 'comp1') object of
#'   'DESeqResults' class, or a dataframe with similar structure. In the case
#'   that a 'DESeqResults' class object is not being used, the object must
#'   contain a 'log2FoldChange' vector/column. Row names must be a subset (but
#'   not necessarily a proper subset) of taxa_names(PS). See Details for more
#'   information.
#' @param comp2 (required) Second comparison (hence 'comp2') object. Refer to
#'   documentation for 'comp1' for remaining details.
#' @param comp1lab Labels for comparison 1. Defaults to c('Decreased Comp1',
#'   'No Change Comp1','Increased Comp1').
#' @param comp2lab Labels for comparison 2. Defaults to c('Decreased Comp2',
#'   'No Change Comp2','Increased Comp2').
#' @param justdata Return only the data table (no plot). Defaults to FALSE.
#' @param PSisRelAbund Does the PS object contain compositional (relative
#'   abundance) taxa counts? Defaults to FALSE.
#'
#' @details
#' \subsection{General}{
#' The results from a single pairwise comparison, such as pre- and
#'   post-treatment, carried out with the DESeq2 package can be plotted or read
#'   and interpreted in tabular form with relative ease. When two pairwise
#'   comparisons are being performed, interpreting the results becomes more
#'   difficult. This function is intended to assist with interpreting the
#'   results from multiple differential abundance analyses carried out with the
#'   DESeq2 r-package. This function takes a phyloseq object ('PS') and two
#'   'DESeqResults' objects ('comp1', 'comp2') and plots the relative abundance
#'   of taxa within 'PS', partitioning the taxa according to their membership
#'   to one the 9 possible cohort combinations determined by their values
#'   specified within 'comp1' and 'comp2'.
#' }
#' \subsection{Approach}{
#' The \code{\link[DESeq2]{DESeq}} function carries out differential abundance
#' testing and produces a 'DESeqDataSet' object. The
#' \code{\link[DESeq2]{results}} function can be used to access the results and
#' create a 'DESeqResults' object, which is a subclass of DataFrame. Note that
#' the 'alpha' parameter for \code{\link[DESeq2]{results}} can be used to
#' specify the significance level of the test being performed. NOTE: Testing
#' with DESeq2 must be carried out, and the non-significant taxa should be
#' removed from the 'comp1' and 'comp2' objects before using this function.
#' Using the log2FoldChange columns in 'comp1' and 'comp2' this function
#' identifies which taxa decrease, do not change, or increase over course of
#' both comparisons. Because there are three options for both comparisons
#' there are 3^2=9 possible combinations, or cohorts, which an OTU may fall
#' into. These cohort assignments are used when plotting the relative abundance
#' plot.
#' }
#'
#'
#' @return A ggplot object.
#'
#' @seealso \code{\link[phyloseq]{plot_bar}}
#'   \code{\link[phyloseq]{transform_sample_counts}}
#'   \code{\link[DESeq2]{DESeq}} \code{\link[DESeq2]{results}}
#'   \code{\link[DESeq2]{DESeqDataSet}}
#'
#' @examples
#' \dontrun{
#' cohort_relabund(
#'   PS=prune_samples(sample_data(WWTP_Impact)$site %in% c(1,2,3,4),
#'                    WWTP_Impact),
#'   comp1=sigtab,
#'   comp2=sigtab.2vs3,
#'   comp1lab=c('Decreased at Effluent',
#'              'No change at Effluent',
#'              'Increased at effluent'),
#'   comp2lab=c('Decreased btwn plants',
#'              'No change btwn plants',
#'              'Increased btwn plants'))
#' }
#'
#' @export

cohort_relabund <- function(PS,xvar='SampleID',taxfill='Phylum',comp1,comp2,comp1lab=c('Decreased Comp1','No Change Comp1','Increased Comp1'),comp2lab=c('Decreased Comp2','No Change Comp2','Increased Comp2'), justdata=FALSE, PSisRelAbund=FALSE){

	ls.comp1.up <- rownames(comp1[comp1$log2FoldChange>0,])
	ls.comp1.down <- rownames(comp1[comp1$log2FoldChange<0,])

	ls.comp2.up <- rownames(comp2[comp2$log2FoldChange>0,])
	ls.comp2.down <- rownames(comp2[comp2$log2FoldChange<0,])

	if (PSisRelAbund==FALSE){
		psra <- transform_sample_counts(PS,function(x) {x / sum(x)})
	} else {
		psra <- PS
	}

	bp <- plot_bar(psra, x=xvar, fill=taxfill)
	bpdat <- bp$data

	#bpdat$comp1 <- comp1lab[2]
	#bpdat$comp1[bpdat$OTU %in% ls.comp1.up]
	bpdat$comp1cat[bpdat$OTU %in% ls.comp1.down] <- comp1lab[1]
	bpdat$comp1cat[!bpdat$OTU %in% c(ls.comp1.down,ls.comp1.up)] <- comp1lab[2]
	bpdat$comp1cat[bpdat$OTU %in% ls.comp1.up] <- comp1lab[3]
	bpdat$comp1cat <-
		factor(bpdat$comp1cat,
					 levels=c(comp1lab[1], comp1lab[2], comp1lab[3]), ordered=TRUE)


	bpdat$comp2cat[bpdat$OTU %in% ls.comp2.down] <- comp2lab[1]
	bpdat$comp2cat[!bpdat$OTU %in% c(ls.comp2.down,ls.comp2.up)] <- comp2lab[2]
	bpdat$comp2cat[bpdat$OTU %in% ls.comp2.up] <- comp2lab[3]
	bpdat$comp2cat <-
		factor(bpdat$comp2cat,
					 levels=c(comp2lab[1], comp2lab[2], comp2lab[3]), ordered=TRUE)

	if (justdata==TRUE){
		bpdat
	} else {

	p1 <-
	ggplot(bpdat, aes_string(x=xvar, y='Abundance', fill=taxfill)) +
		geom_bar(stat='identity', position='stack') +
		theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=0)) +
		facet_grid(comp2cat~comp1cat) +
		theme(legend.position = 'bottom') +
		#theme(legend.title=element_blank()) +
		theme(axis.title.x=element_blank()) #+
		#geom_vline(xintercept=35-5-0.5, linetype='solid', alpha=1.0, size=1.5)
		#ylab('Relative Abundance (Phylum)')

	p1
	}
}


