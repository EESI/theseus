context('quality profile functions')

test_that('qualcontour returns ggplot.',{

  fns <- sort(list.files(system.file('testdata',package='theseus'),full.names=TRUE))
  f_path <- fns[grepl('R1', fns)]
  r_path <- fns[grepl('R2', fns)]

  x <- qualcontour(f_path,r_path,n_samples=2,percentile=.25,nc=1)

  expect_is(x,'ggplot')

})

test_that('prev returns ggplot.',{

  x <- prev(WWTP_Impact)

  expect_is(x,'ggplot')

})