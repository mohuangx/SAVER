
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
              beta = ceiling(post.beta*1000*scale.sf)/1000,
              predicted = ceiling(mu*1000/scale.sf)/1000))
}
