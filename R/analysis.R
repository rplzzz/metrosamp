### Analysis functions for metrosamp results.

#' Analysis functions for MCMC outputs
#'
#' These functions produce summary statistics, such as expectation values or
#' credible intervals, for MCMC objects.
#'
#' @param mcrslt A \code{metrosamp} object, or a list of metrosamp objects
#'   representing multiple chains from the same inference.
#' @param thin Thin mc results to a smaller number of samples before applying the
#' analysis.  This can speed up analyses on large sample runs.
#' @name analysis-functions
NULL

#' @describeIn analysis-functions Maximum A-Posteriori parameters
#' @export
MAP <- function(mcrslt, thin=NULL)
{
    ss <- getsamples(mcrslt, thinto=thin, includelp = TRUE)
    indx <- which.max(ss$lp)
    ss$samples[indx,]
}

#' @describeIn analysis-functions Expectation value for parameters
#' @export
EV <- function(mcrslt, thin=NULL)
{
    samps <- getsamples(mcrslt, thinto=thin)
    apply(samps, 2, mean)
}

#' @describeIn analysis-functions Equal-tailed credibile intervals for parameters
#' @param ci The credible interval to produce (default: 95\%)
#' @export
CI <- function(mcrslt, ci=0.95, thin=NULL)
{
    samps <- getsamples(mcrslt, thinto=thin)
    p <- (1-ci)/2
    q <- c(p, 1-p)
    t(apply(samps, 2, function(x){quantile(x, q)}))
}
