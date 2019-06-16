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


test_that('cor2cov generates the correct covariance matrix', {
    cormat <- matrix(c(1.0, 0.5, -0.5, 0.5, 1.0, 0.0, -0.5, 0.0, 1.0), nrow=3)
    expect_error(cor2cov(cormat, 1.0),'nrow\\(cormat\\) not equal to length\\(scales\\)')

    scales <- c(1.0, 2.0, 3.0)
    covmat <- cor2cov(cormat, scales)

    ## R stats package provides a function to convert covariance to correlation.  Our
    ## function is correct if we can use that to recover our original correlation matrix.
    expect_equal(stats::cov2cor(covmat), cormat)
})

test_that('proposal step distribution is as expected', {
    ngen <- 1000
    scale <- c(1,1)
    pstp <- matrix(nrow=ngen, ncol=2)
    for(i in 1:ngen) {
        pstp[i,] <- proposal_step(scale)
    }
    expect_lt(cor(pstp)[1,2], 0.1)

    sclmat <- matrix(c(1,0.9, 0.9, 1), nrow=2)
    pstpcor <- matrix(nrow=ngen, ncol=2)
    for(i in 1:ngen) {
        pstpcor[i,] <- proposal_step(sclmat)
    }
    expect_gt(cor(pstpcor)[1,2], 0.85)
})
