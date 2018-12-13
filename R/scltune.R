
#' Find a scale vector that produces a desired acceptance rate
#'
#' This function attempts to find a scale vector that will produce a desired
#' acceptance rate for Metropolis sampling of the given probability density
#' funciton.  A secondary effect of this function is that it runs the sampler as
#' it does its work, so the parameter set should be in a high-probability region
#' of the parameter space by the time it is finished.
#'
#' This function is completely experimental, and in its current incarnation
#' represents the most naive way possible of accomplishing its goal.  Feedback
#' on its performance in real-world problems is welcome.
#'
#' @section More details:
#'
#' What we are doing here is using the standard function minimization algorithms
#' of the \code{optim} function to try to find a scale vector that drives the
#' difference between the target acceptance and actual acceptance to zero.  The
#' wrinkle in this strategy is that the actual acceptance rate is a stochastic
#' function of the scale vector, and the minimization algorithms aren't really
#' set up to deal with that.
#'
#' In practice, the algorithms seem to do all right at finding their way to
#' something that is reasonably close to the target value (it probably helps
#' that there are likely many such values; probably they form something like an
#' elliptical surface in the scale factor parameter space).  However, the
#' stochastic behavior seems to mess up the algorithms' convergence criteria
#' pretty badly.  Often the algorithms fly right by a point that achieves the
#' target acceptance \emph{exactly} and settle on one that is rather far off.
#'
#' To combat this tendency, we keep track of the best scale vector we've seen so
#' far (in terms of getting close to the target), and we always report
#' \emph{that} as the result, even if the optimization algorithm actually
#' stopped on something else.  (This is similar to, and inspired by, the
#' way that stochastic gradient descent is used in some machine learning
#' applications).  Ideally we would do this in the optimization algorithm code,
#' but we don't have easy access to that, so instead we keep track of this in
#' the \emph{objective function}, and then fish that information out of the
#' objective function's environment once the optimizer is finished.  It's a
#' little ugly, but it gets the job done.
#'
#' @param lpost Log-posterior function
#' @param p0 Starting parameter values for the sampler.
#' @param scl_guess An initial guess for the scale vector.  If omitted, one will
#' be generated randomly.
#' @param target_accept Target acceptance rate.  Must be between 0 and 1.
#' @param nsamp Number of samples to run each time we evaluate a new scale
#' vector.  Larger values make the estimates of the acceptance probability more
#' robust, but also make the calculation take longer.
#' @return A \code{metrosamp} structure that (hopefully) is tuned for efficient
#' sampling.
#' @export
scltune <- function(lpost, p0, scl_guess=NULL, target_accept=0.25, nsamp=100)
{
    if(is.null(scl_guess)) {
        nvar <- length(p0)
        scl_guess <- runif(nvar, 0.5, 1.5)
    }
    else {
        assertthat::assert_that(all(scl_guess > 0),
                                msg='All components of scl_guess must be > 0.')
    }
    assertthat::assert_that(target_accept > 0 && target_accept < 1)


    ## Create a metropolis sample object
    sampobj <- metrosamp(lpost, p0, 1, 1, scl_guess)

    ## From here on out, we are going to work in terms of the log of the scale
    ## vector.  That allows us to ensure that scale value will be positive,
    ## no matter what the solver algorithm does.
    scl_guess <- log(scl_guess)

    ## Create our objective function
    objfun <- create_objfunc(lpost, sampobj, nsamp, target_accept)

    ## Run the minimization
    optrslt <- optim(scl_guess, objfun, method="BFGS")

    ## Grab the bestx and besty values from the objective function's environment
    bestx <- environment(objfun)[['bestx']]
    besty <- environment(objfun)[['besty']]

    ## If the final acceptance rate is within 0.05 of the target value, call it
    ## good.  (The actual function returns the square of the discrepancy, so
    ## that's what we compare to)
    if(besty > 0.0025) {
        warning('Acceptance rate found by minimization is out of tolerance.  ',
                'bestx= ', bestx, '  besty= ', besty)
    }

    ## We want to retun the sampobj object that we have been carrying around in
    ## the objective function.  We can extract it using `environment`.
    sampobj <- environment(objfun)[['sampobj']]

    ## The code may have kept evaluating after we found our best values, so
    ## reset them to the best values, if necessary
    if(!isTRUE(all.equal(bestx, sampobj$scale))) {
        sampobj <- metrosamp(lpost, sampobj, nsamp, 1, bestx)
    }

    message('optrslt$value = ', optrslt$value, ' accept= ', sampobj$accept,
            ' diff= ', (sampobj$accept-target_accept)^2, ' code= ', optrslt$convergence)

    ## sampobj _should_ have the optimal scale factor already in it, but we'll
    ## make sure of it before we return it.  Don't forget that the parameter
    ## values are logs of the actual scale.
    sampobj$scale <- exp(optrslt$par)

    sampobj
}


## Helper function to create an objective function for scltune.
create_objfunc <- function(lpost, sampobj, nsamp, tgt)
{
    bestx <- NULL
    besty <- Inf
    function(logscl) {
        scl <- exp(logscl)
        ## TODO: should really use the log-post for plast.

        ## This will update the sampobj in our enclosing environment.
        sampobj <<- metrosamp(lpost, sampobj, nsamp, 1, scl)

        message('accept: ', sampobj$accept, '  val: ', (sampobj$accept-tgt)^2)
        y <- (sampobj$accept - tgt)^2
        if(y < besty) {
            bestx <<- scl
            besty <<- y
        }
        y
    }
}
