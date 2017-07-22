#' Calculates SAVER posterior
#'
#' Given prediction and prior variance, calculates the Gamma posterior
#' distribution parameters for a single gene.
#'
#' Let \eqn{\alpha} be the shape parameter and \eqn{\beta} be the rate
#' parameter of the prior Gamma distribution. Then, the posterior Gamma
#' distribution will be
#' \deqn{Gamma(y + \alpha, sf + \beta),}
#' where y is the observed gene count and sf is the size factor.
#'
#' @param y A vector of observed gene counts.
#'
#' @param mu A vector of predictions from \code{\link{expr.predict}}.
#'
#' @param sf Vector of normalized size factors.
#'
#' @param scale.sf Mean of the original size factors.
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression}
#' \item{\code{alpha}}{Posterior Gamma shape parameter}
#' \item{\code{beta}}{Posterior Gamma rate parameter}
#' \item{\code{predicted}}{Regression prediction}
#' \item{\code{lower.95}}{Lower 95\% confidence interval}
#' \item{\code{upper.95}}{Upper 95\% confidence interval}
#'
#' @importFrom stats qgamma
calc.post <- function(y, mu, sf, scale.sf) {
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  if (var(mu) == 0) {
    prior.beta <- rep(calc.b(y, mu, sf)[1], n)
    prior.alpha <- mu*prior.beta
  } else{
    a <- calc.a(y, mu, sf)
    b <- calc.b(y, mu, sf)
    k <- calc.k(y, mu, sf)
    var.method <- which.min(c(a[2], b[2], k[2]))
    if (var.method == 1) {
      prior.alpha <- a[1]
      prior.beta <- a[1]/mu
    } else if (var.method == 2) {
      prior.beta <- b[1]
      prior.alpha <- mu*b[1]
    } else {
      prior.alpha <- mu^2/k[1]
      prior.beta <- mu/k[1]
    }
  }
  post.alpha <- prior.alpha + y
  post.beta <- prior.beta + sf
  lambda.hat <- post.alpha/post.beta
  return(list(estimate = ceiling(lambda.hat*1000/scale.sf)/1000,
              alpha = ceiling(post.alpha*1000)/1000,
              beta = ceiling(post.beta*1000*scale.sf)/1000))
}
