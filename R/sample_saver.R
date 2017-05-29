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
#' data("linnarsson")
#'
#' saver1 <- saver(linnarsson, npred = 5)
#' samp1 <- sample.saver(saver1, seed = 50)

sample.saver <- function(x, rep = 1, efficiency.known = FALSE,
                         seed = NULL) {
  if (!inherits(x, "saver")) {
    stop("x must be a saver object.")
  }
  rep <- as.integer(rep)
  if (rep <= 0) {
    stop("rep must be a positive integer.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (rep == 1) {
    if (efficiency.known) {
      samp <- t(sapply(1:x$ngenes, function(i)
        rnbinom(x$ncells, mu = x$estimate[i, ], size = x$alpha[i, ])))
    } else {
      samp <- t(sapply(1:x$ngenes, function(i)
        rgamma(x$ncells, x$alpha[i, ], x$beta[i, ])))
    }
    rownames(samp) <- x$genes
    colnames(samp) <- x$cells
  } else {
    samp <- vector("list", rep)
    for (j in 1:rep) {
      if (efficiency.known) {
        samp[[j]] <- t(sapply(1:x$ngenes, function(i)
          rnbinom(x$ncells, mu = x$estimate[i, ], size = x$alpha[i, ])))
      } else {
        samp[[j]] <- t(sapply(1:x$ngenes, function(i)
          rgamma(x$ncells, x$alpha[i, ], x$beta[i, ])))
      }
      rownames(samp[[j]]) <- x$genes
      colnames(samp[[j]]) <- x$cells
    }
  }
  return(samp)
}
