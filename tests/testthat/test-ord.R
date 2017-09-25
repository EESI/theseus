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

  X <- data.frame(x=WWTP_Impact@sam_data@.Data[[75]],
                  y=WWTP_Impact@sam_data@.Data[[76]],
                  z=WWTP_Impact@sam_data@.Data[[77]])

  x <- contplot('x','y','z',X,method='linear',removeMissing=FALSE)
  expect_is(x,'ggplot')

  x <- contplot('x','y','z',X,method='spline',removeMissing=FALSE)
  expect_is(x,'ggplot')

  x <- contplot('x','y','z',X,method='loess',removeMissing=FALSE)
  expect_is(x,'ggplot')


})