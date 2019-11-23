#' Optimizes variance
#'
#' Finds the prior parameter that maximizes the marginal likelihood given
#' the prediction.
#'
#' \code{calc.a} returns a prior alpha parameter assuming constant
#' coefficient of variation. \code{calc.b} returns a prior beta parameter
#' assuming constant Fano factor. \code{calc.k} returns a prior variance
#' parameter assuming constant variance.
#'
#' @param y A vector of observed gene counts.
#'
#' @param mu A vector of predictions from \code{\link{expr.predict}}.
#'
#' @param sf Vector of normalized size factors.
#'
#' @return A vector with the optimized parameter and the negative
#' log-likelihood.
#'
#'
#' @importFrom  stats optimize ppoints uniroot var
#'
#' @rdname optimize_variance
#' @export

calc.abk <- function(y, mu, sf, type, samp) {
  n <- length(y)
  if (sum(mu) == 0) {
    return(c(0, 0))
  }
  if (type == "a") {
      interval_max <- var(y/sf)/mean(y/sf)^2
      f <- calc_loglik_a
  } else if (type == "b") {
      interval_max <- var(y/sf)/mean(y/sf)
      f <- calc_loglik_b
  } else {
      interval_max <- var(y/sf)
      f <- calc_loglik_k
  }
  v.vec <- optimize(f, interval = c(0, interval_max), y = y, mu = mu,sf = sf)
  v <- v.vec$minimum
  loglik <- v.vec$objective
  
  min.v <- -f(1e-05, y, mu, sf)
  mle.v <- -f(v, y, mu, sf)
  if (mle.v - min.v < 0.5) {
    max.v <- -f(interval_max, y, mu, sf)
    if (mle.v-max.v > 10) {
      v.max <- uniroot(function(x) f(x, y, mu, sf) + mle.v - 10, c(1e-05, interval_max))$root
    } else {
      v.max <- interval_max
    }
    samp <- samp*v.max
    loglik <- f(samp, y, mu, sf)
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    v <- sum(samp * loglik3)
    loglik <- f(v, y, mu, sf)
  }
  
  
  if (type == "a" || type == "b") {
     return(c(1/v, loglik))
  } else {
     return(c(v, loglik))
  }
  
}