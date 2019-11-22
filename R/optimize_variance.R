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
calc.a <- function(y, mu, sf) {
  n <- length(y)
  if (sum(mu) == 0) {
    return(c(0, 0))
  }
  interval_max <- var(y/sf)/mean(y/sf)^2
  a.vec <- optimize(calc_loglik_a, interval = c(0, interval_max), y = y, mu = mu, sf = sf)
  a <- a.vec$minimum
  a.loglik <- a.vec$objective
  min.a <- -calc_loglik_a(1e-05, y, mu, sf)
  mle.a <- -calc_loglik_a(a, y, mu, sf)
  if (mle.a - min.a < 0.5) {
    max.a <- -calc_loglik_a(interval_max, y, mu, sf)
    if (mle.a-max.a > 10) {
      a.max <- uniroot(function(x) calc_loglik_a(x, y, mu, sf) + mle.a - 10,  c(1e-05, interval_max))$root
    } else {
      a.max <- interval_max
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*a.max
    loglik <- calc_loglik_a(samp, y, mu, sf)
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    a <- sum(samp * loglik3)
    a.loglik <- calc_loglik_a(a, y, mu, sf)
  }
  return(c(1/a, a.loglik))
}

#' @rdname optimize_variance
#' @export
calc.b <- function(y, mu, sf) {
  n <- length(y)
  if (sum(mu) == 0) {
    return(c(0, 0))
  }
  interval_max <- var(y/sf)/mean(y/sf)
  b.vec <- optimize(calc_loglik_b, interval = c(0, interval_max), y = y, mu = mu, sf = sf)
  b <- b.vec$minimum
  b.loglik <- b.vec$objective
  min.b <- -calc_loglik_b(1e-05, y, mu, sf)
  mle.b <- -calc_loglik_b(b, y, mu, sf)
  if (mle.b - min.b < 0.5) {
    max.b <- -calc_loglik_b(interval_max, y, mu, sf)
    if (mle.b-max.b > 10) {
      b.max <- uniroot(function(x) calc_loglik_b(x, y, mu, sf) + mle.b - 10, c(1e-05, interval_max))$root
    } else {
      b.max <- interval_max
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*b.max
    loglik <- calc_loglik_b(samp, y, mu, sf)
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    b <- sum(samp * loglik3)
    b.loglik <- calc_loglik_b(b, y, mu, sf)
  }
  return(c(1/b, b.loglik))
}

#' @rdname optimize_variance
#' @export
calc.k <- function(y, mu, sf) {
  n <- length(y)
  if (sum(mu) == 0) {
    return(c(0, 0))
  }
  interval_max <- var(y/sf)
  k.vec <- optimize(calc_loglik_k, interval = c(0, interval_max), y = y, mu = mu,sf = sf)
  k <- k.vec$minimum
  k.loglik <- k.vec$objective
  min.k <- -calc_loglik_k(1e-05, y, mu, sf)
  mle.k <- -calc_loglik_k(k, y, mu, sf)
  if (mle.k - min.k < 0.5) {
    max.k <- -calc_loglik_k(interval_max, y, mu, sf)
    if (mle.k-max.k > 10) {
      k.max <- uniroot(function(x) calc_loglik_k(x, y, mu, sf) + mle.k - 10, c(1e-05, interval_max))$root
    } else {
      k.max <- interval_max
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*k.max
    loglik <- calc_loglik_k(samp, y, mu, sf)
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    k <- sum(samp * loglik3)
    k.loglik <- calc_loglik_k(k, y, mu, sf)
  }
  return(c(k, k.loglik))
}

#' @rdname optimize_variance
#' @export
calc.abk <- function(y, mu, sf, type) {
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
  interval_max <- var(y/sf)
  k.vec <- optimize(f, interval = c(0, interval_max), y = y, mu = mu,sf = sf)
  k <- k.vec$minimum
  k.loglik <- k.vec$objective
  min.k <- -f(1e-05, y, mu, sf)
  mle.k <- -f(k, y, mu, sf)
  if (mle.k - min.k < 0.5) {
    max.k <- -f(interval_max, y, mu, sf)
    if (mle.k-max.k > 10) {
      k.max <- uniroot(function(x) f(x, y, mu, sf) + mle.k - 10, c(1e-05, interval_max))$root
    } else {
      k.max <- interval_max
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*k.max
    loglik <- f(samp, y, mu, sf)
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    k <- sum(samp * loglik3)
    k.loglik <- f(k, y, mu, sf)
  }
  return(c(k, k.loglik))
}