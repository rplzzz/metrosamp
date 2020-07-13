
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
#' \item{err}{Flag indicating whether an error was caught for the proposal}
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
#' \item{Fix the batch capability, which is currently wrong.}
#' \item{Store the last log-posterior value, so continuation runs don't have to
#' recompute it.}
#' \item{Add code to compute MCSE.}
#' \item{Add option to run functions on MC samples.}
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
#' @param ckpt_name Name of checkpoint files.  Set to \code{NULL} to disable.
#' @param ckpt_freq Frequency to write checkpoint values.  Either a \code{difftime}
#' value, or a numeric value, which will be interpreted as a time in hours. Default
#' is to write a checkpoint every hour.
#' @param lp0 Log-posterior for the starting parameters (p0).  If not supplied
#' it will be calculated automatically.
#' @return A \code{metrosamp} structure of Monte Carlo outputs (described in Details).
#' @importFrom stats runif rnorm
#' @export
metrosamp <- function(lpost, p0, nsamp, batchlen, scale=NULL, debug=FALSE,
                      ckpt_name=NULL, ckpt_freq=1, lp0=NA)
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


    samples <- matrix(nrow=nsamp, ncol=length(p0))
    colnames(samples) <- names(p0)
    prop <- matrix(nrow=nsamp, ncol=length(p0))
    colnames(prop) <- names(p0)
    proplp <- as.numeric(rep(NA, nsamp))
    samplp <- as.numeric(rep(NA, nsamp))
    ratio <- as.numeric(rep(NA, nsamp))
    accept <- rep(0, nsamp)
    errflag <- rep(FALSE, nsamp)

    current_samp <- p0
    if(is.na(lp0)) {
        current_lp <- lpost(p0)
    }
    else {
        current_lp <- lp0
    }

    if(!is.finite(current_lp)) {
        stop('Illegal p0 value:  p0: ', p0, '  log-post: ', current_lp)
    }

    ### Set up checkpointing
    if(is.null(ckpt_name)) {
        do_ckpt <- FALSE
    }
    else {
        do_ckpt <- TRUE
    }

    ## Even if we're not doing checkpoints, we write progress messages at intervals
    ## given by the checkpoint frequency
    if(is.numeric(ckpt_freq)) {
        ckpt_freq <- as.difftime(ckpt_freq, units='hours')
    }
    else if(!inherits(ckpt_freq, 'difftime')) {
        stop('Invalid value for ckpt_freq:  ', ckpt_freq)
    }
    ## record start time
    last_ckpt <- Sys.time()



    ### Ouput function
    metrosamp_struct <- function(i=nsamp) {
        paccept <- sum(accept[1:i])/i
        if(debug) {
            structure(
                list(samples=samples[1:i, ], samplp=samplp[1:i],
                     accept=paccept, plast=current_samp, scale=scale,
                     proposals=prop[1:i, ], proplp=proplp[1:i], ratio=ratio[1:i],
                     prop_accepted=accept[1:i], err=errflag[1:i]),
                class=c('metrosamp', 'list'))
        }
        else {
            structure(
                list(samples=samples[1:i, ], samplp=samplp[1:i],
                     accept=paccept, plast=current_samp, scale=scale),
                class=c('metrosamp','list'))
        }
    }

    for(i in 1:nsamp) {
        if(batchlen == 1) {
            newprop <- current_samp + proposal_step(scale)
            prop[i,] <- newprop
            proplp[i] <- tryCatch(lpost(newprop),error = err)
            if(!is.finite(proplp[i])) {
                if(is.na(proplp[i])) {
                    ## There was an error.  Record it, and set the proposal lp to -Inf
                    errflag[i] <- TRUE
                    proplp[i] <- -Inf
                }
            }
            ratio[i] <- exp(proplp[i] - current_lp)
            if(runif(1) < ratio[i]) {
                ## Accept proposal params
                samples[i,] <- current_samp <- newprop
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
            stop('Batch sampling is not properly implemented.  Do not use')
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

        ## Write progress message and checkpoint if necessary
        tm <- Sys.time()
        runtime <- tm - last_ckpt
        if(runtime >= ckpt_freq) {
            message('time: ', runtime,'  iteration: ', i)
            last_ckpt <- tm
            if(do_ckpt) {
                saveRDS(metrosamp_struct(i), ckpt_name)
            }
        }

    }

    metrosamp_struct()
}

## Error handler
err <- function(e) {
    message(conditionMessage(e))
    NA_real_
}

## proposal generation function
## scale can be either a vector of scale factors or a covariance matrix
proposal_step <- function(scale)
{
    if(is.matrix(scale)) {
        nvar <- nrow(scale)
        MASS::mvrnorm(1, mu=rep(0,nvar), Sigma=scale)
    }
    else {
        nvar <- length(scale)
        scale*rnorm(nvar, 0, 1.0)
    }
}
