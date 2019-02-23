context('Utilities')

set.seed(867-5309)
nocor <- rnorm(100)
allcor <- rep(1, 100)
somecor <- stats::arima.sim(n=100, list(ar=0.5))
samps <- matrix(c(nocor, allcor, somecor), ncol=3)

test_that('Effective size calculation works', {
    ne <- neff(samps)
    expect_equivalent(ne[1], 100)
    expect_equivalent(ne[2], 0)
    ## We don't know exactly what the Neff should be for somecor, since that
    ## depends on random factors, but we can reasonably guess that it should be
    ## somewhere between 30 and 70.  (In fact, for this RNG seed it comes out to
    ## approximately 42.
    expect_gt(ne[3], 30)
    expect_lt(ne[3], 70)
})


