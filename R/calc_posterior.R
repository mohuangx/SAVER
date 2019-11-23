#' Calculates SAVER posterior
#'
#' Given prediction and prior variance, calculates the Gamma posterior
#' distribution parameters for a single gene.
#'x
#' Let \eqn{\alpha} be the shape parameter and \eqn{\beta} be the rate
#' parameter of the prior Gamma distribution. Then, the posterior Gamma
#' distribution will be
#' \deqn{Gamma(y + \alpha, sf + \beta),}
#' where y is the observed gene count and sf is the size factor.
#'
#' @param y A vector of observed gene counts.
#'
#' @param mu A vector of prior means.
#'
#' @param sf Vector of normalized size factors.
#'
#' @param scale.sf Mean of the original size factors.
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression}
#' \item{\code{se}}{Standard error of expression estimate}
#'
#' @importFrom stats qgamma
#' @export
calc.post <- function(y, mu, sf, scale.sf) {
  samp <- (exp((exp(ppoints(100))-1)/2)-1)
  n <- length(y)
  if (length(mu) == 1) {
    mu <- rep(mu, n)
  }
  if (length(sf) == 1) {
    sf <- rep(sf, n)
  }
  if (sum(y) == 0) {
    return(list(estimate = rep(0, n), se = rep(0, n)))
  }
  if (var(mu) == 0) {
    prior.beta <- rep(calc_abk(y, mu, sf, "b", samp)[1], n)
    prior.alpha <- mu*prior.beta
  } else{
    a <- tryCatch(calc_abk(y, mu, sf, "a", samp), error = function(cond) {
      return(c(0, Inf))
    })
    b <- tryCatch(calc_abk(y, mu, sf, "b", samp), error = function(cond) {
      return(c(0, Inf))
    })
    k <- tryCatch(calc_abk(y, mu, sf, "k", samp), error = function(cond) {
      return(c(0, Inf))
    })
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
  se <- sqrt(post.alpha/post.beta^2)
  return(list(estimate = unname(ceiling(lambda.hat*1000*scale.sf)/1000),
              se = unname(ceiling(se*1000*scale.sf)/1000)))
}