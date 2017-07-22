#' Combine function.
#'
#' Combine function for \code{foreach} backend.
#'
#' Genes are split across multiple workers with each worker outputing a list of
#' estimates, alpha, and beta values. This combine function combines the output
#' from the multiple workers into one list.
#'
#' @param LL1 List of SAVER results produced from first worker.
#'
#' @param LL2 List of SAVER results prodcued from another worker.
#'
#' @return A list of combined results


combine.mat <- function(LL1, LL2) {
  estimate <- rbind(LL1$estimate, LL2$estimate)
  alpha <- rbind(LL1$alpha, LL2$alpha)
  beta <- rbind(LL1$beta, LL2$beta)
  return(list(estimate = estimate, alpha = alpha, beta = beta))
}
