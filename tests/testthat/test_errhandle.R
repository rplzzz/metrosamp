context('Handling errors in sampling')

## parameters common to several tests
nvar <- 2
nsamp <- 50
errstep <- 10

## Posterior function that can keep track of the number of times it has been called
## and throw an error every 10th try
counter <- new.env(parent=emptyenv())
counter$i <- 0
errpost <- function(parms) {
    counter$i <- counter$i + 1
    if(counter$i %% errstep == 0) {
        stop('Clear the laundromat!!  This whirl-o-matic just had a nuclear meltdown!!')
    }
    normpost(parms)
}


test_that('Errors are caught and reported.', {
    out <- expect_message(metrosamp(errpost, rep(1, nvar), nsamp, 1, rep(1,nvar), debug=TRUE),
                          regexp='laundromat')
    ## Every 10th call is an error, but we do one extra call before sampling starts
    ## to get the initial log-posterior.
    expect_equal(out$err, seq(2,nsamp+1) %% errstep == 0)

    ## Check that we still get messages, even when debug is false
    out <- expect_message(metrosamp(errpost, rep(1, nvar), nsamp, 1, rep(1,nvar), debug=FALSE),
                          regexp='laundromat')

})


