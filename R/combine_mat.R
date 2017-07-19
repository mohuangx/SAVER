


combine.mat <- function(LL1, LL2) {
  estimate <- rbind(LL1$estimate, LL2$estimate)
  alpha <- rbind(LL1$alpha, LL2$alpha)
  beta <- rbind(LL1$beta, LL2$beta)
  return(list(estimate = estimate, alpha = alpha, beta = beta))
}
