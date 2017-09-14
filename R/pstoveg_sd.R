#' @import phyloseq
NULL

#' converts the otu_table slot of a phyloseq object to a vegan-compatible matrix
#'
#' physeq2veg_otu is a helper function intended to convert the species/taxa count slot of a phyloseq object to a vegan-friendly matrix. 
#'
#' @param PS (required) a phyloseq object
#'
#' @return A matrix containing a phyloseq object's otu_table slot. 
#'
#' @seealso \code{\link[phyloseq]{phyloseq-class}} \code{\link[phyloseq]{otu_table-class}} \code{\link[phyloseq]{otu_table}}
#'
#' @examples
#' \dontrun{
#' library(phyloseq)
#' data("GlobalPatterns")
#' # inspect otu_table()
#' dim(otu_table(GlobalPatterns))
#' str(otu_table(GlobalPatterns))
#' taxa_are_rows(GlobalPatterns)
#' gp.otu <- physeq2veg_otu(GlobalPatterns)
#' dim(gp.otu)
#' str(gp.otu)
#' }
#'
#' @export

pstoveg_otu <- function(PS) {
  OTU <- otu_table(PS)
  if (taxa_are_rows(OTU)) {
    OTU <- t(OTU)
  }
  return(as(OTU, "matrix"))
}