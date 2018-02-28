#' Calculate maximum correlation
#'
#' Calculates the maximum absolute correlation between two matrices along the
#' columns
#'
#' This function calculates the maximum absolute correlation for each column of
#' \code{x2} against each column of \code{x1}. The matrices are named and if
#' the names overlap between \code{x1} and \code{x2}, the correlation between
#' the same named entries is set to zero.
#'
#' @param x1 Named matrix 1
#' @param x2 Named matrix 2
#'
#' @return A vector of maximum absolute correlations
#'
#' @examples
#' x1 <- matrix(rnorm(500), 100, 5)
#' x2 <- x1 + matrix(rnorm(500), 100, 5)
#' colnames(x1) <- c("A", "B", "C", "D", "E")
#' colnames(x2) <- c("A", "B", "C", "D", "E")
#' cor(x1, x2)
#' calc.maxcor(x1, x2)
#'
#' @export

calc.maxcor <- function(x1, x2) {
  cormat <- suppressWarnings(cor(x1, x2))
  cormat[is.na(cormat)] <- 0
  if (nrow(cormat) < ncol(cormat)) {
    cind <- sapply(rownames(cormat), function(x)
      which(x == colnames(cormat)))
    rind <- which(sapply(cind, length) == 1)
    cind <- unlist(cind)
  } else {
    rind <- sapply(colnames(cormat), function(x)
      which(x == rownames(cormat)))
    cind <- which(sapply(rind, length) == 1)
    rind <- unlist(rind)
  }
  cormat[cbind(rind, cind)] <- 0
  apply(cormat, 2, function(x) max(abs(x)))
}