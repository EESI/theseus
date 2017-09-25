context('helper functions')

test_that('get_top_taxa returns taxa names of correct length',{

  top <- get_top_taxa(WWTP_Impact,tax_level='Phylum',tax_n=7)

  expect_length(top,7)
  expect_is(top,'character')
  expect_true(all(!is.na(top)))

})


test_that('rename_column_to_other returns taxa names of correct length',{

  df <- phyloseq::psmelt(WWTP_Impact)
  top_taxa <- get_top_taxa(WWTP_Impact)
  other <- table(rename_column_to_other(df,tax_level='Phylum',top_taxa)$Phylum,useNA='always')

  expect_length(other,9)
  expect_equal(unname(other[which(is.na(names(other)))]),0)

})

test_that('pstoveg_otu returns matrix',{

  otu <- pstoveg_otu(WWTP_Impact)

  expect_is(otu,'matrix')

})