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

    ## Test neff for metrosamp objects.
    ngen <- 100
    p0 <- c(1,1)
    scale <- c(0.25, 0.25)

    ## generate metrosamp runs with and without debugging
    set.seed(867-5309)
    ms1 <- metrosamp(rosenbrock, p0, ngen, 1, scale)
    ms2 <- metrosamp(rosenbrock, p0, ngen, 1, scale)
    expect_equal(neff(ms1), neff(ms1$samples))
    expect_equal(neff(list(ms1,ms2)), neff(rbind(ms1$samples, ms2$samples)))
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

test_that('sample extraction works', {
    mk_test_ms <- function(vals, samplp, paccept=0, scale=numeric(0)) {
        if(is.matrix(vals)) {
            samples <- vals
        }
        else {
            samples <- matrix(1, nrow=length(vals), ncol=3) * vals
        }
        if(length(scale)==0) {
            scale <- rep(1, ncol(samples))
        }
        structure(
            list(samples=samples, samplp=samplp, accept=paccept, plast=samples[nrow(samples),], scale=scale),
            class=c('metrosamp','list'))
    }

    ## test single mc struct
    mcruns1 <- mk_test_ms(c(1,2,3), c(-1,-2,-3), 1)

    samps1 <- getsamples(mcruns1)
    expected_samps <- matrix(rep(c(1,2,3), 3), ncol=3)
    expect_is(samps1, 'matrix')
    expect_equal(samps1, expected_samps)

    samps1lp <- getsamples(mcruns1,includelp=TRUE)
    expect_is(samps1lp, 'list')
    expect_equal(samps1lp, list(samples=expected_samps, lp=c(-1,-2,-3)))

    ## Passing a matrix or a list of matrices in should fail
    expect_error(getsamples(mcruns1$samples), 'must be a metrosamp object')
    expect_error(getsamples(list(mcruns1$samples, mcruns1$samples)), 'must be a metrosamp object')

    ### test multiple mc struct.  Give each one 5 copies of the same sample
    mcruns3 <- mapply(mk_test_ms, rep(c(1,2,3), rep(5,3)), rep(c(-1,-2,-3), rep(5,3)), SIMPLIFY=FALSE)
    samps3 <- getsamples(mcruns3)
    expected_samps <- matrix(rep(rep(c(1,2,3), rep(5,3)), 3), ncol=3)
    expect_is(samps3, 'matrix')
    expect_equal(samps3, expected_samps)

    samps3lp <- getsamples(mcruns3, includelp=TRUE)
    expect_is(samps3lp, 'list')
    expect_equal(samps3lp, list(samples=expected_samps, lp=rep(c(-1,-2,-3), rep(5,3))))

    ### test with thinning in effect
    samps3thinlp <- getsamples(mcruns3, 3, TRUE)
    expected_samps <- matrix(rep(c(1,2,3), 3), ncol=3)   # should get one from each chain
    expect_is(samps3thinlp, 'list')
    expect_equal(samps3thinlp, list(samples=expected_samps, lp=c(-1,-2,-3)))

    ### different thinning factors
    samps3thin5 <- getsamples(mcruns3, 5, TRUE)
    expect_equal(samps3thin5,
                 list(samples=samps3[c(3,6,9,12,15),], lp=samps3lp$lp[c(3,6,9,12,15)]))
    samps3thin7 <- getsamples(mcruns3, 6)
    expect_equal(samps3thin7, samps3[c(2,4,6,8,10,12),])
})

test_that('Lists of metrosamp objects can be converted to coda objects', {
    ngen <- 100
    p0 <- c(1,1)
    scale <- c(0.25, 0.25)
    set.seed(867-5309)

    ms1 <- metrosamp(rosenbrock, p0, ngen, 1, scale)
    ms2 <- metrosamp(rosenbrock, p0, ngen, 1, scale)

    c1 <- metrosamp2coda(ms1)
    c12 <- metrosamp2coda(list(ms1,ms2))

    expect_s3_class(c1, 'mcmc')
    expect_s3_class(c12, 'mcmc.list')

    expect_equal(coda::niter(c1), ngen)
    expect_equal(coda::niter(c12), ngen)

    expect_equal(coda::nvar(c1), 2)
    expect_equal(coda::nvar(c12), 2)

    expect_equal(coda::nchain(c1), 1)
    expect_equal(coda::nchain(c12), 2)

    c1t <- metrosamp2coda(ms1, 25)
    expect_s3_class(c1t, 'mcmc')
    expect_equal(coda::niter(c1t), 25)

    c12t <- metrosamp2coda(list(ms1,ms2), 25)
    expect_s3_class(c12t, 'mcmc.list')
    expect_equal(coda::niter(c12t), 25)
    expect_equal(coda::nchain(c12t), 2)
})
