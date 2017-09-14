#' Calculates marginal likelihood
#'
#' Calculates the marginal likelihood given the prediction under constant
#' coefficient of variation (a), Fano factor (b), and variance (k).
#'
#' \code{calc.loglik.a} returns the shifted negative log-likelihood under
#' constant coefficient of variation.
#' \code{calc.loglik.b} returns the shifted negative log-likelihood under
#' constant Fano factor.
#' \code{calc.loglik.k} returns the shifted negative log-likelihood under
#' constant variance.
#'
#' @param a,b,k Prior parameter.
#'
#' @param y A vector of observed gene counts.
#'
#' @param mu A vector of predictions from \code{\link{expr.predict}}.
#'
#' @param sf Vector of normalized size factors.
#'
#' @return A shifted negative marginal log-likelihood.
#'
#'
#' @rdname calc_loglik
#' @export
calc.loglik.a <- function(a, y, mu, sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  func1 <- n/a*log(1/a)
  func2 <- -sum(1/a*log(mu))
  func3 <- -n*lgamma(1/a)
  func4 <- sum(lgamma(y+1/a))
  func5 <- -sum((y+1/a)*log(sf+1/(a*mu)))
  return(-sum(func1, func2, func3, func4, func5))
}

#' @rdname calc_loglik
#' @export
calc.loglik.b <- function(b, y, mu, sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  func1 <- sum(mu/b*log(1/b))
  func2 <- -sum(lgamma(mu/b))
  func3 <- sum(lgamma(y+mu/b))
  func4 <- -sum((y+mu/b)*log(sf+1/b))
  return(-sum(func1, func2, func3, func4))
}

#' @rdname calc_loglik
#' @export
calc.loglik.k <- function(k, y, mu, sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  func3 <- sum(mu^2*log(mu)/k)
  func4 <- -sum(mu^2*log(k)/k)
  func5 <- -sum(lgamma(mu^2/k))
  func6 <- sum(lgamma(y+mu^2/k))
  func7 <- -sum((y+mu^2/k)*log(sf+mu/k))
  return(-sum(func3, func4, func5, func6, func7))
}
