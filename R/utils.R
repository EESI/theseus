#' @importFrom magrittr %>%
#' @export
magrittr::`%>%`

#' @importFrom phyloseq tax_table<-
#' @export
phyloseq::`tax_table<-`

#' @importFrom phyloseq otu_table<-
#' @export
phyloseq::`otu_table<-`

#' @importFrom phyloseq sample_data<-
#' @export
phyloseq::`sample_data<-`

if(getRversion() >= "2.15.1")  utils::globalVariables(c("."))