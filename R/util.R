#### Utility functions

#' Correct Markov chain sample count for autocorrelation.
#'
#' The output of a Markov chain sampling process is autocorrelated; therefore,
#' the effective number of samples is less than the actual number of samples
#' collected.  This function computes the corrected value.
#'
#' The input can be any of:
#' \itemize{
#' \item{A matrix of samples, with samples in rows and variables in columns}
#' \item{A metrosamp object}
#' \item{A list of metrosamp objects}
#' }
#'
#' This function is currently just a thin wrapper around
#' \code{\link[coda]{effectiveSize}}.  As noted below, the method used in this
#' calculation is not ideal, so eventually we will replace this with something
#' more sophisticated.
#'
#'
#' @section Note:
#'
#' The effective sample size is computed as \deqn{N_e =
#' \frac{N}{1+2\sum_{k=1}^\infty} \rho(k),} where \eqn{\rho(k)} is the
#' autocorrelation at lag \eqn{k}.  We calculate this independently for each
#' sampled variable and return the result as a vector.  Generally, if these
#' values differ greatly, you will want the smallest value reported.
#'
#' Andrew Gelman makes the following observation about this calculation:
#'
#' \dQuote{[This definition] for effective sample size doesn’t quite work
#' in practice because (a) you can’t sum to infinity, and (b) it will be too
#' optimistic for chains that haven’t mixed. We have an effective sample size
#' estimate that addresses both these concerns. It’s in formulas (11.7) and
#' (11.8) in Section 11.5 of BDA3.
#'
#' If you don’t have a copy of BDA, it’s also in the Stan manual, on pages
#' 353-354 of the manual for Stan version 2.14.0.}
#' (\url{https://www.johndcook.com/blog/2017/06/27/effective-sample-size-for-mcmc/#comment-933045})
#'
#' We have gone with the simpler definition for now, since this package is
#' mostly intended for light-duty work, but we should consider putting in the
#' more sophisticated version at some point.
#'
#' @param mcmcrslt MCMC results.  See details.
#' @return Vector of N_e values, one for each variable
#' @export
neff <- function(mcmcrslt)
{
    if(is.matrix(mcmcrslt)) {
        samps <- mcmcrslt
    }
    else {
        ## Assume we're dealing with a metrosamp object or a list of same
        samps <- getsamples(mcmcrslt)
    }
    coda::effectiveSize(samps)
}
#' Convert a correlaiton matrix into a covariance matrix
#'
#' Given a correlation matrix an a vector of standard deviations for the individual
#' variables, produce a covariance matrix.
#'
#' This is a convenience function, meant for producing covariance matrices for
#' proposal distributions in cases where you have an idea of what the correlation
#' should be, and you want to specify the scale independently.  No checks are
#' performed to see, for example, whether the correlation matrix is valid.
#'
#' @param cormat Correlation matrix
#' @param scales Vector of scale factors (i.e., standard deviations)
#' @export
cor2cov <- function(cormat, scales)
{
    assertthat::assert_that(is.matrix(cormat))
    assertthat::assert_that(is.vector(scales))
    assertthat::assert_that(nrow(cormat) == ncol(cormat))
    assertthat::assert_that(nrow(cormat) == length(scales))

    ## We do this with two calls to sweep.  One gets the rows, the other
    ## gets the columns
    out <- sweep(cormat, 1, scales, "*", check.margin=FALSE)
    sweep(out, 2, scales, "*", check.margin=FALSE)
}


#' Concatenate metrosamp objects from continuation runs into a single object
#'
#' This function takes \code{\link{metrosamp}} objects from an initial run and
#' metrosamp objects from a continuation run and combines them into a single
#' metrosamp object for the entire collection of iterations.  The arguments can
#' be either two individual metrosamp objects, or two lists of metrosamp
#' objects.  In the latter case, the objects are concatenated pairwise to
#' produce a list of objects for the concatenated runs.
#'
#' Debuggging information from the runs is kept if it is present in both
#' objects; otherwise it is dropped.
#'
#' The intent of this function is to concatenate runs with subsequent runs that
#' started where the first run left off, but no check is made to ensure that the
#' second run really is a continuation of the first.  Concatenating two runs
#' that aren't continuations may result in an error, or in bogus results.
#'
#' @param mslist1 First metrosamp object or list of metrosamp objects to concatenate.
#' @param mslist2 Second metrosamp object or list of metrosamp objects to concatenate.
#' @return A metrosamp object, or a list of metrosamp objects, representing the concatenated runs.
#' @export
concat <- function(mslist1, mslist2)
{
    if(inherits(mslist1, 'metrosamp') && inherits(mslist2, 'metrosamp')) {
        ## Got individual metrosamp objects instead of mslists.  Return the
        ## two objects, concatenated
        concat_single(mslist1, mslist2)
    }
    else {
        ## Concatenate term by term
        mapply(concat_single, mslist1, mslist2, SIMPLIFY = FALSE)
    }
}

## Helper function for concatenate function.
concat_single <- function(ms1, ms2)
{
    assertthat::assert_that(inherits(ms1, 'metrosamp'))
    assertthat::assert_that(inherits(ms2, 'metrosamp'))

    nsamp1 <- nrow(ms1$samples)
    nsamp2 <- nrow(ms2$samples)
    nsamptot <- nsamp1 + nsamp2

    ## Deal with the items that are always present
    samps <- rbind(ms1$samples, ms2$samples)
    samplp <- c(ms1$samplp, ms2$samplp)
    paccept <- (nsamp1*ms1$accept + nsamp2*ms2$accept) / nsamptot
    plast <- ms2$plast
    scale <- ms2$scale

    ## Check to see if debugging is on in _both_ objects.  Debugging information
    ## only makes sense if it was turned on all the way through, so drop it if it's
    ## missing in _either_ object.
    if('proposals' %in% names(ms1) && 'proposals' %in% names(ms2)) {
        props <- rbind(ms1$proposals, ms2$proposals)
        proplp <- c(ms1$proplp, ms2$proplp)
        ratio <- c(ms1$ratio, ms2$ratio)
        prop_accept <- c(ms1$prop_accepted, ms2$prop_accepted)
        structure(
            list(samples=samps, samplp=samplp, accept=paccept, plast=plast, scale=scale,
                 proposals=props, proplp=proplp, ratio=ratio, prop_accepted=prop_accept),
            class=c('metrosamp', 'list'))
    }
    else {
        structure(
            list(samples=samps, samplp=samplp, accept=paccept, plast=plast, scale=scale),
            class=c('metrosamp', 'list'))
    }
}

#' Extract matrix of samples from a \code{metrosamp} object
#'
#' Given a metrosamp object, or a list of \code{metrosamp} objects, extract the
#' matrix of parameter samples.  Optionally, also extract the log-posterior
#' values.
#'
#' If the input is a list of \code{metrosamp} objects, then their samples are
#' merged into a single grand matrix of results.  Doing this discards the
#' information about which samples came from which Markov chains.
#'
#' @param mcruns A \code{metrosamp} object or a list of \code{metrosamp} objects.
#' @param thinto Total number of samples to thin the sample array to.  Many analysis
#' operations scale as \eqn{O(N)} or \eqn{O(N \ln N)}, so thinning can speed these
#' up pretty substantially, especially for runs with poor sampling efficiency.
#' @param includelp If true, return the log-posterior values in a list with the sample
#' matrix.
#' @return If \code{includelp} is \code{FALSE}, a matrix of samples, otherwise, a list
#' with the sample matrix and the vector of log-posterior values.
#' @export
getsamples <- function(mcruns, thinto=NULL, includelp=FALSE)
{
    if(inherits(mcruns, 'metrosamp')) {
        samps <- mcruns$samples
        if(includelp) {
            lp <- mcruns$samplp
        }
    }
    else {
        ## presumably a list of metrosamp objects
        if(!is.list(mcruns) || ! all(sapply(mcruns, inherits, what='metrosamp'))) {
            stop('mcruns argument must be a metrosamp object or list of metrosamp objects')
        }
        samps <- do.call(rbind, lapply(mcruns, function(mco){mco$samples}))
        if(includelp) {
            lp <- do.call(c, lapply(mcruns, function(mco){mco$samplp}))
        }
    }

    if(is.numeric(thinto) && thinto > 0) {
        ntot <- nrow(samps)
        fac <- floor(ntot / thinto)
        keep <- seq(1,ntot) %% fac == 0
        samps <- samps[keep,][1:thinto,]  # second indexing is for cases where thinto does not evenly divide ntot
        if(includelp) {
            lp <- lp[keep][1:thinto]
        }
    }

    if(includelp) {
        list(samples=samps, lp=lp)
    }
    else {
        samps
    }
}

#' Convert a list of metrosamp structures to a \code{\link[coda]{mcmc.list}}
#' structure.
#'
#' This conversion allows all of the analysis and diagnostic functions from the
#' coda package to be applied to metrosamp results
#'
#' @param mslist A list of metrosamp structures
#' @param size Size of the mcmc objects in output coda structure. The results will be
#' thinned as necessary to get to this size.
#' @export
metrosamp2coda <- function(mslist, size=NA) {


    if(inherits(mslist, 'metrosamp')) {
        ## we got only a single metrosamp structure, so return a coda::mcmc object
        if(is.na(size)) {
            thin <- 1
        }
        else {
            thin <- ceiling(nrow(mslist$samples) / size)
        }

        samps <- mslist$samples
        idx <- seq(1,nrow(samps))
        samps <- samps[idx %% thin == 0, ]
        coda::mcmc(samps, thin=thin)
    }
    else {
        if(is.na(size)) {
            thin <- 1
        }
        else {
            thin <- ceiling(nrow(mslist[[1]]$samples) / size)
        }

        coda::mcmc.list(
            lapply(mslist, function(ms) {
                samps <- ms$samples
                idx <- seq(1,nrow(samps))
                samps <- samps[idx %% thin == 0, ]
                coda::mcmc(samps, thin=thin)
            })
        )
    }
}
