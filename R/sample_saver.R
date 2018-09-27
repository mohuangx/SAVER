#' Samples from SAVER
#'
#' Samples from the posterior distribution output by SAVER.
#'
#' The SAVER method outputs a posterior distribution, which we can sample from
#' for downstream analysis. The posterior distribution accounts for uncertainty
#' in the SAVER estimation procedure. If the efficiency is known,
#' negative binomial sampling is performed; otherwise, gamma sampling is
#' performed.
#'
#' @param x A \code{saver} object.
#'
#' @param rep Number of sampled datasets. Default is 1.
#'
#' @param efficiency.known Whether the efficiency is known. Default is
#' \code{FALSE}.
#'
#' @param seed seed used in \code{set.seed}.
#'
#' @return A matrix of expression values sampled from SAVER posterior. If
#' \code{rep} > 1, a list of matrices is returned
#'
#' @examples
#' data("linnarsson_saver")
#'
#' samp1 <- sample.saver(linnarsson_saver, seed = 50)
#' 
#'
#' @importFrom stats rgamma rnbinom
#' @export

sample.saver <- function(x, rep = 1, efficiency.known = FALSE,
                         seed = NULL) {
  if (!is.null(x$alpha)) {
    return(sample.saver.old(x, rep, efficiency.known, seed))
  }
  ncells <- ncol(x$estimate)
  ngenes <- nrow(x$estimate)
  cell.names <- colnames(x$estimate)
  gene.names <- rownames(x$estimate)
  rep <- as.integer(rep)
  if (rep <= 0) {
    stop("rep must be a positive integer.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (rep == 1) {
    if (efficiency.known) {
      samp <- t(sapply(1:ngenes, function(i)
        rnbinom(ncells, mu = x$estimate[i, ],
                size = x$estimate[i, ]^2/x$se[i, ]^2)))
    } else {
      samp <- t(sapply(1:ngenes, function(i)
        rgamma(ncells, x$estimate[i, ]^2/x$se[i, ]^2,
               x$estimate[i, ]/x$se[i, ]^2)))
      samp <- round(samp, 3)
    }
    rownames(samp) <- gene.names
    colnames(samp) <- cell.names
  } else {
    samp <- vector("list", rep)
    for (j in 1:rep) {
      if (efficiency.known) {
        samp[[j]] <- t(sapply(1:ngenes, function(i)
          rnbinom(ncells, mu = x$estimate[i, ],
                  size = x$estimate[i, ]^2/x$se[i, ]^2)))
      } else {
        samp[[j]] <- t(sapply(1:ngenes, function(i)
          rgamma(ncells, x$estimate[i, ]^2/x$se[i, ]^2,
                 x$estimate[i, ]/x$se[i, ]^2)))
        samp[[j]] <- round(samp[[j]], 3)
      }
      rownames(samp[[j]]) <- gene.names
      colnames(samp[[j]]) <- cell.names
    }
  }
  return(samp)
}
