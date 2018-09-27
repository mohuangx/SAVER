#' Calculates gene-to-gene and cell-to-cell SAVER correlation
#'
#' Adjusts for SAVER estimation uncertainty by calculating and adjusting
#' gene-to-gene and cell-to-cell correlation matrices
#'
#' The SAVER estimates that are produced have varying levels of uncertainty
#' depending on the gene and the cell. These functions adjust the gene-to-gene
#' and cell-to-cell correlations of the SAVER estimates to reflect the
#' estimation uncertainty.
#'
#' @param x A \code{saver} object.
#'
#' @param cor.mat If a correlation matrix of the SAVER estimates was already
#' obtained, then it can be provided as an input to avoid recomputation.
#'
#' @return An adjusted correlation matrix.
#'
#' @examples
#' data("linnarsson_saver")
#' 
#' gene.cor <- cor.genes(linnarsson_saver)
#'
#' @importFrom stats cor
#' @export
#' @rdname cor_adjust
cor.genes <- function(x, cor.mat = NULL) {
  if (!is.null(x$alpha)) {
    return(cor.genes.old(x, cor.mat))
  }
  if (is.null(cor.mat)) {
    message("Calculating correlation matrix...")
    cor.mat <- cor(t(x$estimate))
  }
  ngenes <- nrow(x$estimate)
  adj.vec <- rep(0, ngenes)
  for (i in 1:ngenes) {
    adj.vec[i] <- sqrt(var(x$estimate[i, ], na.rm = TRUE)/
                         (var(x$estimate[i, ], na.rm = TRUE) +
                            mean(x$se[i, ]^2, na.rm = TRUE)))
  }
  adj.mat <- outer(adj.vec, adj.vec)
  cor.adj <- adj.mat*cor.mat
  return(cor.adj)
}

#' @export
#' @rdname cor_adjust
cor.cells <- function(x, cor.mat = NULL) {
  if (!is.null(x$alpha)) {
    return(cor.cells.old(x, cor.mat))
  }
  if (is.null(cor.mat)) {
    message("Calculating correlation matrix...")
    cor.mat <- cor(x$estimate)
  }
  ncells <- ncol(x$estimate)
  adj.vec <- rep(0, ncells)
  for (i in 1:ncells) {
    adj.vec[i] <- sqrt(var(x$estimate[, i], na.rm = TRUE)/
                         (var(x$estimate[, i], na.rm = TRUE) +
                            mean(x$se[, i]^2, na.rm = TRUE)))
  }
  adj.mat <- outer(adj.vec, adj.vec)
  cor.adj <- adj.mat*cor.mat
  return(cor.adj)
}

