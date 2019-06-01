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


test_that('Debug information is added when requested.',
      {
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
})
