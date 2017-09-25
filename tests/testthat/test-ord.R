context('ordination functions')

test_that('envtoverlay returns ggplot for PCA and CA.',{

  x <- envtoverlay(WWTP_Impact, c('Mg','Na','K'), ordmet='PCA')
  expect_is(x,'ggplot')

  x <- envtoverlay(WWTP_Impact, c('Mg','Na','K'), ordmet='CA')
  expect_is(x,'ggplot')

})

test_that('constord returns ggplot for various params.',{

  x <- constord(WWTP_Impact,~Mg + Na,method='RDA',scaling=2,tax_level='Phylum',tax_n=7)
  expect_is(x,'ggplot')

  x <- constord(WWTP_Impact,~Mg + Na,method='CCA',scaling=2,tax_level='Phylum',tax_n=7)
  expect_is(x,'ggplot')

  x <- constord(WWTP_Impact,~Mg + Na,method='CCA',scaling=1,tax_level='Phylum',tax_n=7)
  expect_is(x,'ggplot')

  x <- constord(WWTP_Impact,~Mg + Na,method='CCA',facets=~day,scaling=1,tax_level='Phylum',tax_n=7)
  expect_is(x,'ggplot')

})

test_that('constord returns ggplot for various params.',{

  X <- data.frame(phyloseq::sample_data(WWTP_Impact),stringsAsFactors=FALSE)

  x <- contplot('log_Si','log_Cu','log_Mn',X,method='linear',removeMissing=FALSE)
  expect_is(x,'ggplot')

  x <- contplot('log_Si','log_Cu','log_Mn',X,method='spline',removeMissing=FALSE)
  expect_is(x,'ggplot')

  x <- contplot('log_Si','log_Cu','log_Mn',X,method='loess',removeMissing=FALSE)
  expect_is(x,'ggplot')


})