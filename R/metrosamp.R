
#' Perform Metropolis sampling with batch averaging
#'
#' This function performs a basic Metropolis sampling of a user-supplied
#' log-posterior function.  The sampling is done in batches, with the batch
#' means returned as the output.  Setting the batch length to 1 will produce
#' unbatched samples.
#'
#' The output \code{metrosamp} structure will be a list with the following elements:
#' \describe{
#' \item{samples}{Matrix (nsamp x nparam) of parameter samples}
#' \item{samplp}{Vector (nsamp) of log-posterior values for the samples}
#' \item{accept}{Probability of accepting a proposal, averaged across all
#' samples}
#' \item{plast}{Last parameter set.  Can be used to continue sampling where the
#' last run left off.}
#' \item{scale}{Scale factor used in the calculation.  Also useful for
#' continuing a run.}
#' }
#' If the \code{debug} flag is set, the output will have some additional
#' elements that can be used to diagnose the sampling procedure.  If batch
#' sampling is in use, then most of these will pertain to the \emph{last}
#' proposal evaluated in each batch.  Intermediate proposals within a batch are
#' not returned.  Therefore, when debugging proposals, it is best to use a batch
#' length of 1.
#' \describe{
#' \item{proposals}{The proposal parameters evaluated by the sampler.}
#' \item{proplp}{Log-posterior for the proposals.  The same notes apply as to
#' the \code{proposals} entry.}
#' \item{prop_accepted}{Flag indicating whether each proposal was accepted.}
#' }
#'
#' A run can be continued by passing the \code{metrosamp} structure from the
#' previous run as the \code{p0} argument.  If this is done, then the
#' \code{scale} parameter may be omitted, and the new run will use the same
#' scale as the old.  If a scale parameter \emph{is} supplied, then it will
#' override the scale parameter stored in the old structure.
#'
#' @section To Do:
#'
#' \itemize{
#' \item{Add code to compute MCSE.}
#' \item{Add code to compute Neff.}
#' \item{Add option to run functions on MC samples.}
#' \item{Allow covariance matrix for scale parameter.}
#' }
#'
#' @param lpost Log-posterior function
#' @param p0 Starting parameters for sampling, OR a \code{metrosamp} structure
#' from a previous run.
#' @param nsamp Number of batches to run
#' @param batchlen Number of samples per batch
#' @param scale MC step scaler; this will be multiplied by a vector of standard
#' normal deviates to get the proposal step.  Optional if a \code{metrosamp}
#' structure was supplied for \code{p0}; required otherwise.
#' @param debug Flag to turn on additional debugging information.
#' @param lp0 Log-posterior for the starting parameters (p0).  If not supplied
#'   it will be calculated automatically.
#' @return A \code{metrosamp} structure of Monte Carlo outputs (described in Details).
#' @export
metrosamp <- function(lpost, p0, nsamp, batchlen, scale=NULL, debug=FALSE, lp0=NA)
{

    if(inherits(p0, 'metrosamp')) {
        if(is.null(scale)) {
            ## only use the old run's scale if no scale was supplied
            scale <- p0$scale
        }
        p0 <- p0$plast
    }

    assertthat::assert_that(!is.null(scale),
                            msg='If p0 is not a metrosamp object, then scale must be supplied.')

    nvar <- length(p0)

    samples <- matrix(nrow=nsamp, ncol=length(p0))
    prop <- matrix(nrow=nsamp, ncol=length(p0))
    proplp <- as.numeric(rep(NA, nsamp))
    samplp <- as.numeric(rep(NA, nsamp))
    ratio <- as.numeric(rep(NA, nsamp))
    accept <- rep(0, nsamp)

    current_samp <- p0
    if(is.na(lp0)) {
        current_lp <- lpost(p0)
    }
    else {
        current_lp <- lp0
    }
    for(i in 1:nsamp) {
        if(batchlen == 1) {
            prop[i,] <- current_samp + scale*rnorm(nvar, 0, 1.0)
            proplp[i] <- lpost(prop[i,])
            ratio[i] <- exp(proplp[i] - current_lp)
            if(runif(1) < ratio[i]) {
                ## Accept proposal params
                samples[i,] <- current_samp <- prop[i,]
                samplp[i] <- current_lp <- proplp[i]
                accept[i] <- 1
            }
            else {
                ## reject proposal params
                samples[i,] <- current_samp
                samplp[i] <- current_lp
            }
        }
        else {
            ## sample a batch
            message('batch: ', i)
            batch <- metrosamp(lpost, current_samp, batchlen, 1, scale, debug, current_lp)
            ## set this iteration's result using the result of the batch
            samples[i,] <- batchmean <- apply(batch$samples, 2, mean)
            samplp[i] <- lpost(batchmean)     # give us the actual log posterior for the batch average params
            accept[i] <- mean(batch$accept)
            ## just set the proposal to whatever the last proposal was
            if(debug) {
                prop[i,] <- batch$proposals[batchlen,]
                proplp[i] <- batch$proplp[batchlen]
                ratio[i] <- batch$ratio[batchlen]
            }
        }
    }
    paccept <- sum(accept)/nsamp
    if(debug) {
        structure(
            list(samples=samples, proposals=prop, proplp=proplp, samplp=samplp, ratio=ratio,
                 prop_accepted=accept, accept=paccept, plast=samples[nsamp,], scale=scale),
            class=c('metrosamp', 'list'))
    }
    else {
        structure(
            list(samples=samples, samplp=samplp, accept=paccept, plast=samples[nsamp,], scale=scale),
            class=c('metrosamp','list'))
    }
}
