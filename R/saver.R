#' Runs SAVER
#'
#' Recovers expression using the SAVER method.
#'
#' The SAVER method starts by estimating the prior mean and variance for the
#' true expression level for each gene and cell. The prior mean is obtained
#' through predictions from a LASSO Poisson regression for each gene
#' implemented using the \code{glmnet} package. Then, the variance is estimated
#' through maximum likelihood assuming constant variance, Fano factor, or
#' coefficient of variation variance structure for each gene. The posterior
#' distribution is calculated and the posterior mean is reported as the SAVER
#' estimate.
#'
#' @param x An expression count matrix. The rows correspond to genes and
#' the columns correspond to cells.
#'
#' @param size.factor Vector of cell size normalization factors.
#' If \code{x} is already normalized or normalization is not desired, use
#' \code{size.factor = 1}. Default uses mean library size normalization.
#'
#' @param npred Number of genes for regression prediction. Selects the top
#' \code{npred} genes in terms of mean expression for regression prediction.
#' Default is all genes.
#'
#' @param pred.genes Indices of specific genes to perform regression
#' prediction. Overrides \code{npred}. Default is all genes.
#'
#' @param pred.genes.only Return expression levels of only \code{pred.genes}.
#' Default is FALSE (returns expression levels of all genes).
#'
#' @param parallel If \code{TRUE}, uses the parallel backend \code{foreach} to
#' parallelize prediction. Must register parallel backend beforehand
#' (\code{doParallel}, \code{doMC}). See example below.
#'
#' @param dfmax Maximum number of genes used in prediction model. Default is
#' 300.
#'
#' @param nfolds Number of folds to be used in cross validation. Default is 5.
#'
#' @param remove.zero.genes Whether to remove genes with all zero. Default is FALSE.
#'
#' @param cutoff If the model with the lowest mean cross-validation error has
#' Poisson deviance greater than the mean cross-validation error of the null
#' model subtracted by the estimated standard error of the mean
#' cross-validation error of the null model times \code{cutoff}, then the null
#' model is selected. Default is 0.
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression}
#' \item{\code{alpha}}{Posterior Gamma shape parameter}
#' \item{\code{beta}}{Posterior Gamma rate parameter}
#' \item{\code{predicted}}{Regression prediction}
#' \item{\code{lower.95}}{Lower 95\% confidence interval}
#' \item{\code{upper.95}}{Upper 95\% confidence interval}
#' \item{\code{size.factor}}{Size factors used to normalize expression.}
#' \item{\code{ngenes}}{Number of genes.}
#' \item{\code{ncells}}{Number of cells.}
#' \item{\code{genes}}{Gene names.}
#' \item{\code{cells}}{Cell names.}
#'
#' @examples
#' data("linnarsson")
#'
#' # predictions for top 5 highly expressed genes
#' saver1 <- saver(linnarsson, npred = 5)
#'
#' # predictions for certain genes
#' genes <- c("Thy1", "Mbp", "Stim2", "Psmc6", "Rps19")
#' genes.ind <- which(rownames(linnarsson) %in% genes)
#' saver2 <- saver(linnarsson, pred.genes = genes.ind)
#'
#' # return only certain genes
#' saver3 <- saver(linnarsson, pred.genes = genes.ind, pred.genes.only = TRUE)
#'
#' # Parallel
#' \dontrun{
#' require(doParallel)
#' registerDoParallel(cores = 2)
#' system.time(saver(linnarsson, npred = 20))
#' system.time(saver(linnarsson, npred = 20, parallel = TRUE))
#' }
#'
#' @export
saver <- function(x, size.factor = NULL, npred = NULL, pred.genes = NULL,
                  pred.genes.only = FALSE, parallel = FALSE, dfmax = 300,
                  nfolds = 5, remove.zero.genes = FALSE, cutoff = 0) {
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  if (min(colSums(x)) == 0) {
    x <- x[, colSums(x) != 0]
    message("Removing cells with zero expression.")
  }
  nonzero <- which(rowSums(x) != 0)
  if (remove.zero.genes) {
    x <- x[nonzero, ]
    message("Removing genes with zero expression.")
  }
  np <- dim(x)
  ngenes <- as.integer(np[1])
  ncells <- as.integer(np[2])
  if (is.null(size.factor)) {
    sf <- colSums(x)/mean(colSums(x))
    scale.sf <- 1
  } else if (length(size.factor) == ncells) {
    sf <- size.factor/mean(size.factor)
    scale.sf <- mean(size.factor)
  } else if (size.factor == 1) {
    sf <- rep(1, ncells)
    scale.sf <- 1
  } else if (min(size.factor) <= 0) {
    stop("Size factor must be greater than 0")
  } else {
    stop("Not a valid size factor")
  }
  if (!is.null(pred.genes)) {
    if (!is.integer(pred.genes) | min(pred.genes) < 1 |
        max(pred.genes) > ngenes) {
      stop("pred.genes must be row indices of x")
    }
    npred <- length(pred.genes)
  } else if (is.null(npred)) {
    npred <- length(nonzero)
    pred.genes <- nonzero
  } else if (npred < ngenes) {
    pred.genes <- order(rowMeans(x), decreasing = TRUE)[1:npred]
  } else {
    stop("npred must be less than number of rows in x")
  }
  if (length(nonzero) < ngenes) {
    pred.genes <- pred.genes[pred.genes %in% nonzero]
    x.norm <- sweep(x[nonzero, ], 2, sf, "/")
  } else {
    x.norm <- sweep(x, 2, sf, "/")
  }
  x.est <- log(sweep(x.norm, 2, 1/sf, "+"))
  npred <- length(pred.genes)
  if (pred.genes.only)
    ngenes <- npred
  gene.means <- rowMeans(x.norm)
  message("calculating predictions...")
  mu <- matrix(0, 5, ncells)
  if (npred > 5) {
    s <- sample(1:npred, 5)
    t1 <- Sys.time()
    for (i in 1:5) {
      mu[i, ] <- expr.predict(t(x.est[-pred.genes[s[i]], ]),
                              x[pred.genes[s[i]], ]/sf, dfmax, nfolds)
    }
    t2 <- Sys.time()
    for (i in 1:5) {
      post <- calc.post(x[pred.genes[s[i]], ], mu[i, ], sf, scale.sf)
    }
    t3 <- Sys.time()
    t.diff1 <- (t2-t1)/5
    t.diff2 <- (t3-t2)/5
    units(t.diff1) <- "secs"
    units(t.diff2) <- "secs"
  }
  if (parallel & nworkers > 1) {
    if (npred > 5) {
      nworkers <- foreach::getDoParWorkers()
      t3 <- t.diff1*npred/nworkers*1.1 + t.diff2*ngenes/nworkers*1.1
      units(t3) <- "mins"
      message("Approximate finish time: ", t2+t3)
    }
    message("Running in parallel: ", nworkers, " workers")
    gene.list <- chunk2(pred.genes, nworkers)
    mu.par <- foreach::foreach(i = 1:nworkers, .combine = rbind,
                           .packages = c("glmnet", "SAVER")) %dopar% {
      mu.temp <- matrix(0, length(gene.list[[i]]), ncells)
      for (j in 1:length(gene.list[[i]])) {
        mu.temp[j, ] <- expr.predict(t(x.est[-gene.list[[i]][j], ]),
                                     x[gene.list[[i]][j], ]/sf, dfmax,
                                     nfolds, seed = gene.list[[i]][j],
                                     cutoff)
      }
    return(mu.temp)
    }
  } else {
    if (parallel & nworkers == 1) {
      message("Only one worker assigned! Running sequentially...")
    }
    if (npred > 5) {
      t3 <- t.diff1*npred*1.1 + t.diff2*ngenes*1.1
      units(t3) <- "mins"
      message("Approximate finish time: ", t2+t3)
    }
    mu.par <- matrix(0, npred, ncells)
    for (i in 1:length(pred.genes)) {
      mu.par[i, ] <- expr.predict(t(x.est[-pred.genes[i], ]),
                                  x[pred.genes[i], ]/sf, dfmax, nfolds,
                                  seed = pred.genes[i], cutoff)
    }
  }
  message("Predictions finished. Calculating posterior...")
  if (pred.genes.only) {
    out <- lapply(1:6, function(x) matrix(0, npred, ncells))
    mu <- mu.par
    for (i in 1:length(pred.genes)) {
      post <- calc.post(x[pred.genes[i], ], mu[i, ], sf, scale.sf)
      out[[1]][i, ] <- post$estimate
      out[[2]][i, ] <- post$alpha
      out[[3]][i, ] <- post$beta
      out[[5]][i, ] <- post$lower.95
      out[[6]][i, ] <- post$upper.95
    }
    out[[4]] <- mu
    gene.names <- rownames(x)[pred.genes]
    cell.names <- colnames(x)
  } else {
    out <- lapply(1:6, function(x) matrix(0, ngenes, ncells))
    mu <- matrix(gene.means, ngenes, ncells)
    mu[pred.genes, ] <- mu.par
    for (i in 1:ngenes) {
      post <- calc.post(x[i, ], mu[i, ], sf, scale.sf)
      out[[1]][i, ] <- post$estimate
      out[[2]][i, ] <- post$alpha
      out[[3]][i, ] <- post$beta
      out[[5]][i, ] <- post$lower.95
      out[[6]][i, ] <- post$upper.95
    }
    out[[4]] <- mu
    gene.names <- rownames(x)
    cell.names <- colnames(x)
  }
  out.named <- lapply(out, function(x) {rownames(x) <- gene.names;
                                        colnames(x) <- cell.names; x})
  out.named[[7]] <- sf*scale.sf
  out.named[[8]] <- ngenes
  out.named[[9]] <- ncells
  out.named[[10]] <- gene.names
  out.named[[11]] <- cell.names
  names(out.named[[7]]) <- cell.names
  names(out.named) <- c("estimate", "alpha", "beta", "predicted", "lower.95",
                        "upper.95", "size.factor", "ngenes", "ncells",
                        "genes", "cells")
  class(out.named) <- "saver"
  message("Done!")
  return(out.named)
}


