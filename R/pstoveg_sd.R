#' converts the sam_data slot of a phyloseq object to a vegan-compatible matrix
#'
#' physeq2veg_sd is a helper function intended to convert the sample data slot
#' of a phyloseq object to a vegan-friendly matrix.
#'
#' @param PS (required) a phyloseq object
#'
#' @return A matrix containing a phyloseq object's sam_data slot.
#'
#' @seealso \code{\link[phyloseq]{phyloseq-class}}
#'   \code{\link[phyloseq]{sample_data-class}}
#'   \code{\link[phyloseq]{sample_data}}
#'
#' @examples
#' \dontrun{
#' library(theseus)
#' library(phyloseq)
#' data(WWTP_Impact, package='theseus')
#' dim(sample_data(WWTP_Impact))
#' sampdat <- pstoveg_sd(WWTP_Impact)
#' dim(sampdat)
#'
#' data(GlobalPatterns, package='phyloseq')
#' dim(sample_data(GlobalPatterns))
#' sampdat.gp <-pstoveg_sd(GlobalPatterns)
#' dim(sampdat.gp)
#'
#' # move altered sample data back to phyloseq
#' sampdat.altered <- sampdat
#' sampdat.altered$TotDisP_PercentMax <- vegan::decostand(sampdat$TotDisP,
#'                                                        method='max')
#' sample_data(wwtp) <- as.data.frame(sampdat.altered)
#' }
#'
#' @export

pstoveg_sd <- function(PS) {
  sd <- sample_data(PS)
  return(as(sd,"data.frame"))
}
