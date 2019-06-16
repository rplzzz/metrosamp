context('Basic sampling functions')

## parameters common to several tests
nvar <- 2
nsamp <- 1000

test_that('Unbatched sampling works.', {
    out <- metrosamp(normpost, rep(1, nvar) , 1000, 1, rep(1, nvar))

    expect_s3_class(out,'metrosamp')
    expect_type(out,'list')

    expect_is(out$samples, 'matrix')
    expect_equal(dim(out$samples), c(nsamp, nvar))

    expect_true(is.vector(out$samplp))
    expect_equal(length(out$samplp), nsamp)

    expect_equivalent(out$samples[nsamp,], out$plast)

    ## Check debug info is NOT included
    expect_false('proposals' %in% names(out))
    expect_false('proplp' %in% names(out))
    expect_false('ratio' %in% names(out))
    expect_false('prop_accepted' %in% names(out))

})


test_that('Debug information is added when requested.', {
    out <- metrosamp(normpost, rep(1, nvar), 1000, 1, rep(1,nvar), debug=TRUE)

    expect_s3_class(out,'metrosamp')
    expect_type(out,'list')

    expect_is(out$proposals, 'matrix')
    expect_equal(dim(out$proposals), c(nsamp, nvar))

    expect_true(is.vector(out$proplp))
    expect_equal(length(out$proplp), nsamp)

    expect_true(is.vector(out$ratio))
    expect_equal(length(out$ratio), nsamp)

    expect_true(is.vector(out$prop_accepted))
    expect_equal(length(out$prop_accepted), nsamp)

})


test_that('Batched sampling works.', {
    out <- metrosamp(normpost, rep(1, nvar) , 1000, 10, rep(1, nvar))

    expect_s3_class(out,'metrosamp')
    expect_type(out,'list')

    expect_is(out$samples, 'matrix')
    expect_equal(dim(out$samples), c(nsamp, nvar))

    expect_true(is.vector(out$samplp))
    expect_equal(length(out$samplp), nsamp)
})


test_that('Increasing the proposal scale decreases acceptance probability.', {
    set.seed(8675309)
    scl <- rep(1, nvar)
    out1 <- metrosamp(normpost, rep(1, nvar), 1000, 1, scl)

    set.seed(8675309)
    out2 <- metrosamp(normpost, rep(1, nvar), 1000, 1, 4*scl)

    expect_lt(out2$accept, out1$accept)
})

test_that('Correlation matrix for proposals works.', {
    set.seed(867-5309)
    ## run an uncorrelated version for comparison
    scl_uncor <- rep(1, nvar)
    out1 <- metrosamp(normpost, rep(1,nvar), 1000, 1, scl_uncor, debug=TRUE)
    expect_lt(cor(out1$proposals)[1,2], 0.1)   ## proposal correlation should be small

    scl_cor <- diag(nrow=nvar, ncol=nvar)
    scl_cor[1,2] <- scl_cor[2,1] <- 0.9
    out2 <- metrosamp(normpost, rep(1, nvar), 1000, 1, scl_cor, debug=TRUE)
    ## The proposals entry in the output gives us the actual proposals, not the
    ## proposal steps.  Because the proposals are influenced by the posterior, which
    ## is symmetric, they will not be as correlated as the step distribution.  Still,
    ## they should be noticeably more correlated than the uncorrelated version.
    expect_gt(cor(out2$proposals)[1,2], 0.7)
    ## The correlated version should also accept more often.  This is a bit counterintuitive,
    ## but the important thing here is that the correlation doesn't increase the scale, along
    ## the major axis (relative to the uncorrelated version), but it _does_ decrease it
    ## along the minor axis.  Therefore, the samples must be more concentrated near the center.
    ## Given our symmetric posterior, this will lead to higher acceptance rates.
    expect_gt(out2$accept, out1$accept)

    ## Now check that scaling the correlated step size larger does decrease the acceptance rate.
    sclfac <- rep(2,nvar)
    sclfac[1] <- 3   # make it asymmetric
    scl_cor2 <- cor2cov(scl_cor, sclfac)
    out3 <- metrosamp(normpost, rep(1, nvar), 1000, 1, scl_cor2)
    expect_lt(out3$accept, out2$accept)
})


test_that('Continuation runs work.', {
    set.seed(8675309)
    out1 <- metrosamp(normpost, rep(1, nvar) , 1000, 1, rep(1, nvar))

    ## continue with same scale
    set.seed(8675309)
    out2 <- metrosamp(normpost, out1, 1000, 1)
    ## Check tht we get valid output
    expect_s3_class(out1,'metrosamp')
    expect_is(out1$samples, 'matrix')
    expect_equal(dim(out1$samples), c(nsamp, nvar))

    expect_true(is.vector(out1$samplp))
    expect_equal(length(out1$samplp), nsamp)

    expect_equal(out1$scale, out2$scale)

    ## Try overriding the scale
    set.seed(8675309)
    out3 <- metrosamp(normpost, out2, 1000, 1, rep(0.5, nvar))
    expect_gt(out3$accept, out1$accept)
    expect_gt(out3$accept, out2$accept)

    ## Check that omitting the scale produces an error if it's not a
    ## continuation run.
    expect_error(metrosamp(normpost, rep(1, nvar), 1000, 1))
})


test_that('Sampling a simple posterior produces the expected distribution.', {
    set.seed(8675309)
    out <- metrosamp(normpost, rep(1, nvar), 1000, 10, rep(1,nvar))

    distmean <- apply(out$samples, 2, mean)
    distsd <- apply(out$samples, 2, sd)

    ## As a crude estimate, we _should_ be able to get the mean of the
    ## distribution right to within 1/sqrt(N)
    tol <- 1/sqrt(nsamp)
    expect_equal(distmean, rep(1, nvar), tolerance=tol) # expected mean is 1,1

    ## The higher the moment we are looking at, the less accurate it's going to
    ## be for a given number of samples.  In test calculations, the standard
    ## deviations of the samples _almost_ come in at 1/sqrt(N) accuracy.  Give
    ## them a little bit of slack for these automated tests.
    tol <- tol*1.25
    expect_equal(distsd, rep(0.25, nvar), tolerance=tol)
})

test_that('Sampler preserves attributes', {
    vn <- c('a','b','c','d')
    lpost <- function(p) {
        stopifnot(!is.null(names(p)) && all(names(p) == vn))
        stopifnot(!is.null(attr(p, 'bogon flux')) && attr(p, 'bogon flux') == 0.0)
        -sum((p - c(1.0,2.0,3.0,4.0))^2)
    }

    p0 <- c(1.0,2.0,3.0,4.0)
    sig <- rep(1,4)
    expect_error({
        metrosamp(lpost, p0, 100, 1, sig)
    })
    names(p0) <- vn
    expect_error({
        metrosamp(lpost, p0, 100, 1, sig)
    })
    attr(p0, 'bogon flux') <- 2.5   # Nonzero bogon flux!
    expect_error({
        metrosamp(lpost, p0, 100, 1, sig)
    })
    attr(p0, 'bogon flux') <- 0.0   # Bogon flux nominal
    expect_silent({
        ms <- metrosamp(lpost, p0, 100, 1, sig)
    })
    expect_silent({
        ms2 <- metrosamp(lpost, ms, 100, 1)
    })

    ## try it again with a correlated proposal step
    sigmat <- diag(nrow=4, ncol=4)
    sigmat[1,2] <- sigmat[2,1] <- 0.9
    attr(p0, 'bogon flux') <- 2.5
    expect_error({
        metrosamp(lpost, p0, 100, 1, sigmat)
    })
    attr(p0, 'bogon flux') <- 0.0
    expect_silent({
        ms <- metrosamp(lpost, p0, 100, 1, sigmat)
    })
    expect_silent({
        ms2 <- metrosamp(lpost, ms, 100, 1)
    })
})
