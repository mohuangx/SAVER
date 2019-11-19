
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
  func5 <- -sum(gammaln(mu^2/k))
  func6 <- sum(gammaln(y+mu^2/k))
  func7 <- -sum((y+mu^2/k)*log(sf+mu/k))
  return(-sum(func3, func4, func5, func6, func7))
}
