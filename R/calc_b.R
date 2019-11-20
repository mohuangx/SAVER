#' @rdname optimize_variance
#' @export
calc.b <- function(y, mu, sf) {
  n <- length(y)
  if (sum(mu) == 0) {
    return(c(0, 0))
  }
  var_over_mean <- var(y/sf)/mean(y/sf)
  b.vec <- optimize(calc_loglik_b, interval = c(0, var_over_mean), y = y, mu = mu, sf = sf)
  b <- b.vec$minimum
  b.loglik <- b.vec$objective
  min.b <- -calc_loglik_b(1e-05, y, mu, sf)
  mle.b <- -calc_loglik_b(b, y, mu, sf)
  max.b <- -calc_loglik_b(var_over_mean, y, mu, sf)
  if (mle.b - min.b < 0.5) {
    if (mle.b-max.b > 10) {
      b.max <- uniroot(function(x) calc_loglik_b(x, y, mu, sf) + mle.b - 10, c(1e-05, var_over_mean))$root
    } else {
      b.max <- var_over_mean
    }
    samp <- (exp((exp(ppoints(100))-1)/2)-1)*b.max
    loglik <- mapply(calc_loglik_b, b = samp, MoreArgs = list(y = y, mu = mu, sf = sf))
    loglik2 <- exp(-loglik-min(-loglik))
    loglik3 <- loglik2/sum(loglik2)
    b <- mean(sample(samp, 10000, replace = TRUE, prob = loglik3))
    b.loglik <- calc_loglik_b(b, y, mu, sf)
  }
  return(c(1/b, b.loglik))
}