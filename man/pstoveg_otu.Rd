% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pstoveg_otu.R
\name{pstoveg_otu}
\alias{pstoveg_otu}
\title{converts the otu_table slot of a phyloseq object to a vegan-compatible matrix}
\usage{
pstoveg_otu(PS)
}
\arguments{
\item{PS}{(required) a phyloseq object}
}
\value{
A matrix containing a phyloseq object's otu_table slot.
}
\description{
physeq2veg_otu is a helper function intended to convert the species/taxa
count slot of a phyloseq object to a vegan-friendly matrix. This function
ensures that sites/samples are rows and species are columns.
}
\examples{
\dontrun{
library(theseus)
library(phyloseq)
data(WWTP_Impact, package='theseus')
dim(otu_table(WWTP_Impact))
taxa_are_rows(WWTP_Impact)
otu <- pstoveg_otu(WWTP_Impact)
dim(otu)

data(GlobalPatterns, package='phyloseq')
dim(otu_table(GlobalPatterns))
taxa_are_rows(GlobalPatterns)
otu.gp <-pstoveg_otu(GlobalPatterns)
dim(otu.gp)

# move transformed OTU table back to phyloseq
wwtp <- WWTP_Impact
otu.ra <- vegan::decostand(otu, method='total')
otu_table(wwtp) <- otu_table(otu.ra,
                             taxa_are_rows = taxa_are_rows(WWTP_Impact))
}

}
\seealso{
\code{\link[phyloseq]{phyloseq-class}}
  \code{\link[phyloseq]{otu_table-class}}
  \code{\link[phyloseq]{otu_table}}
}
