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

test_that('concatenating metrosamp objects works', {
    ngen <- 100
    p0 <- c(1,1)
    scale <- c(0.25, 0.25)

    ## generate metrosamp runs with and without debugging
    set.seed(867-5309)
    ms1 <- metrosamp(rosenbrock, p0, ngen, 1, scale, FALSE)
    ms2 <- metrosamp(rosenbrock, ms1, ngen, 1, debug=FALSE)
    set.seed(867-5309)
    ms1d <- metrosamp(rosenbrock, p0, ngen, 1, scale, TRUE)
    ms2d <- metrosamp(rosenbrock, ms1d, ngen, 1, debug=TRUE)

    ## now generate the equivalent single runs
    set.seed(867-5309)
    ms12 <- metrosamp(rosenbrock, p0, 2*ngen, 1, scale, FALSE)
    set.seed(867-5309)
    ms12d <- metrosamp(rosenbrock, p0, 2*ngen, 1, scale, TRUE)

    ## concatenate the first group of runs
    ms12c <- concat(ms1, ms2)
    ms12cd <- concat(ms1d, ms2d)
    ms12mix <- concat(ms1d, ms2)
    ms12mix2 <- concat(ms1, ms2d)

    expect_equal(ms12c, ms12)
    expect_equal(ms12cd, ms12d)
    expect_equal(ms12mix, ms12)
    expect_equal(ms12mix2, ms12)
})

test_that('concatenating lists of metrosamp objects works', {
    ngen <- 100
    p0 <- c(1,1)
    scale <- c(0.25, 0.25)

    set.seed(867-5309)
    ms1a <- metrosamp(rosenbrock, p0, ngen, 1, scale)
    ms1b <- metrosamp(rosenbrock, ms1a, ngen, 1)
    ms2a <- metrosamp(rosenbrock, p0, ngen, 1, scale)
    ms2b <- metrosamp(rosenbrock, ms2a, ngen, 1)

    set.seed(867-5309)
    ms1 <- metrosamp(rosenbrock, p0, 2*ngen, 1, scale)
    ms2 <- metrosamp(rosenbrock, p0, 2*ngen, 1, scale)

    mslc <- concat(list(ms1a, ms2a), list(ms1b, ms2b))
    msl <- list(ms1, ms2)

    expect_equal(mslc, msl)
})
