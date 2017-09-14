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
#' @rdname optimize_variance
#' @export
calc.a <- function(y, mu, sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  a.vec <- optimize(calc.loglik.a, interval = c(0, var(y/sf)/mean(y/sf)^2),
                    y = y, mu = mu, sf = sf)
  a <- a.vec$minimum
  a.loglik <- a.vec$objective
  min.a <- -calc.loglik.a(1e-05, y, mu, sf)
  mle.a <- -calc.loglik.a(a, y, mu, sf)
  max.a <- -calc.loglik.a(var(y/sf)/mean(y/sf)^2, y, mu, sf)
  if (mle.a - min.a < 0.5) {
    if (mle.a-max.a > 10) {
      a.max <- uniroot(function(x) calc.loglik.a(x, y, mu, sf) + mle.a - 10,
                       c(1e-05, var(y/sf)/mean(y/sf)^2))$root
    } else {
      a.max <- var(y/sf)/mean(y/sf)^2
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*a.max
    loglik <- mapply(calc.loglik.a, a = samp,
                     MoreArgs = list(y = y,
                                     mu = mu,
                                     sf = sf))
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    a <- mean(sample(samp, 10000, replace = TRUE, prob = loglik3))
    a.loglik <- calc.loglik.a(a, y, mu, sf)
  }
  return(c(1/a, a.loglik))
}

#' @rdname optimize_variance
#' @export
calc.b <- function(y, mu, sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  if (sum(mu) == 0) {
    return(c(0, 0))
  }
  b.vec <- optimize(calc.loglik.b, interval = c(0, var(y/sf)/mean(y/sf)), y = y, mu = mu,
                    sf = sf)
  b <- b.vec$minimum
  b.loglik <- b.vec$objective
  min.b <- -calc.loglik.b(1e-05, y, mu, sf)
  mle.b <- -calc.loglik.b(b, y, mu, sf)
  max.b <- -calc.loglik.b(var(y/sf)/mean(y/sf), y, mu, sf)
  if (mle.b - min.b < 0.5) {
    if (mle.b-max.b > 10) {
      b.max <- uniroot(function(x) calc.loglik.b(x, y, mu, sf) + mle.b - 10,
                       c(1e-05, var(y/sf)/mean(y/sf)))$root
    } else {
      b.max <- var(y/sf)/mean(y/sf)
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*b.max
    loglik <- mapply(calc.loglik.b, b = samp,
                     MoreArgs = list(y = y,
                                     mu = mu,
                                     sf = sf))
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    b <- mean(sample(samp, 10000, replace = TRUE, prob = loglik3))
    b.loglik <- calc.loglik.b(b, y, mu, sf)
  }
  return(c(1/b, b.loglik))
}

#' @rdname optimize_variance
#' @export
calc.k <- function(y, mu, sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  k.vec <- optimize(calc.loglik.k, interval = c(0, var(y/sf)), y = y, mu = mu,
                    sf = sf)
  k <- k.vec$minimum
  k.loglik <- k.vec$objective
  min.k <- -calc.loglik.k(1e-05, y, mu, sf)
  mle.k <- -calc.loglik.k(k, y, mu, sf)
  max.k <- -calc.loglik.k(var(y/sf), y, mu, sf)
  if (mle.k - min.k < 0.5) {
    if (mle.k-max.k > 10) {
      k.max <- uniroot(function(x) calc.loglik.k(x, y, mu, sf) + mle.k - 10,
                       c(1e-05, var(y/sf)))$root
    } else {
      k.max <- var(y/sf)
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*k.max
    loglik <- mapply(calc.loglik.k, k = samp,
                     MoreArgs = list(y = y,
                                     mu = mu,
                                     sf = sf))
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    k <- mean(sample(samp, 10000, replace = TRUE, prob = loglik3))
    k.loglik <- calc.loglik.k(k, y, mu, sf)
  }
  return(c(k, k.loglik))
}
