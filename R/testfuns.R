#### Test functions for exercising the sampler

#' Calculate the Rosenbrock function
#'
#' The Rosenbrock function is a function that is constructed to be difficult to
#' minimize due to its shape.  In two dimensions the function is shaped like a
#' long, narrow, curving valley.  The global minimum is at the lowest point in
#' the valley.  It is the dramatic difference in slope between the directions
#' along and across the valley that makes it so hard to find the minimum.
#' Therefore, it is often used as a test of optimizers or of Monte Carlo
#' Samplers.
#'
#' In two dimensions the definition of the Rosenbrock function is
#' \deqn{(1-x)^2 + 100(y-x^2)^2.}
#' In higher dimensions this generalizes to
#' \deqn{\sum_{i=1}^{N-1} 100(x_{i+1} - x_i^2)^2 + (1-x_i)^2.}
#' The 2-D case has a minimum at \eqn{(1,1)}.  The 3-D case has a similar
#' minimum at \eqn{(1,1,1)}.  For \eqn{4\le N\le 7}, there is still a global
#' minimum at \eqn{(1,1, \ldots, 1)}, but there is also a second local minimum
#' at \eqn{(-1, 1, 1, \ldots, 1)}.
#'
#' @param x Vector of input values.  The dimension of the function will be
#' inferred from the length of this vector.
#' @export
rosenbrock <- function(x)
{
    N <- length(x)
    stopifnot(N >= 2)

    y <- x[2:N]
    x <- x[1:(N-1)]

    sum( (1-x)^2 + 100*(y-x^2)^2 )
}



#' @describeIn rosenbrock Log-transformed Rosenbrock function, suitable for use as a log-posterior.
#' @export
logrosen <- function(x)
{
    -log1p(rosenbrock(x))
}


#' Log-posterior for a multi-variate normal distribution
#'
#' This is just about the simplest test you can run on a sampler.  The
#' distribution is centered at \eqn{(1,1,\ldots,1)}
#'
#' @param x Vector of input values.
#' @export
normpost <- function(x)
{
    -2*sum((x-1)^2)
}
