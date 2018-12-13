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
