get_top_taxa <- function(PS,tax_level='Phylum',tax_n=7){

  OTU <- phyloseq::otu_table(PS)@.Data %>%
    data.frame(stringsAsFactors = FALSE)

  if (PS@otu_table@taxa_are_rows) OTU <- t(OTU)

  TAX <- phyloseq::tax_table(PS)@.Data %>%
    data.frame(stringsAsFactors = FALSE) %>%
    dplyr::mutate_(id=~rownames(.))

  top <- data.frame(abundance=colSums(OTU),id=colnames(OTU),stringsAsFactors=FALSE) %>%
    dplyr::left_join(TAX,by='id') %>%
    dplyr::group_by_(tax_level) %>%
    dplyr::summarize_(abundance=~sum(abundance)) %>%
    dplyr::arrange_(~dplyr::desc(abundance)) %>%
    as.data.frame()

  top <- top[1:tax_n,1]

  return(top)

}

rename_column_to_other <- function(df,tax_level,top_taxa){

  df[[tax_level]] <- as.character.factor(df[[tax_level]])
  df[[tax_level]][!(df[[tax_level]] %in% top_taxa)] <- 'Other'
  df[[tax_level]] <- factor(df[[tax_level]],levels=c(top_taxa,'Other'))

  return(df)

}

get_quant <- function(x,y,q) x[which(cumsum(y)/sum(y) >= q)][[1]]

qual_scores <- function(path,percentile=.25,forward=TRUE){

  srqa <- ShortRead::qa(path)

  df <- srqa[['perCycle']]$quality

  quant <- as.vector(by(df,df$Cycle,function(x) get_quant(x$Score,x$Count,percentile),simplify=TRUE))

  if (forward){
    out <- matrix(rep(quant,length(quant)),nrow=length(quant))
  }else{
    out <- matrix(rep(quant,length(quant)),nrow=length(quant),byrow=TRUE)
  }

  return(out)

}