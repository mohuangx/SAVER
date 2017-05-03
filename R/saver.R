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
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered normalized expression}
#' \item{\code{alpha}}{Posterior Gamma shape parameter}
#' \item{\code{beta}}{Posterior Gamma rate parameter}
#' \item{\code{predicted}}{Regression prediction}
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
#' require(doParallel)
#' registerDoParallel(cores = 2)
#' system.time(saver(linnarsson, npred = 20))
#' system.time(saver(linnarsson, npred = 20, parallel = TRUE))
#'
#' @export


saver <- function(x, size.factor = NULL, npred = NULL, pred.genes = NULL,
                  pred.genes.only = FALSE, parallel = FALSE, dfmax = 300,
                  nfolds = 5) {
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  if (min(rowSums(x)) == 0 | min(colSums(x)) == 0)
    stop("Please remove genes or cells with all zero expression")
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
        max(pred.genes) > as.integer(np[1])) {
      stop("pred.genes must be row indices of x")
    }
    npred <- length(pred.genes)
  } else if (is.null(npred)) {
    npred <- ngenes
    pred.genes <- 1:ngenes
  } else if (npred < ngenes) {
    pred.genes <- order(rowMeans(x), decreasing = TRUE)[1:npred]
  } else {
    stop("npred must be less than number of rows in x")
  }
  if (pred.genes.only)
    ngenes <- npred
  x.norm <- sweep(x, 2, sf, "/")
  x.est <- log(sweep(x.norm, 2, 1/sf, "+"))
  gene.means <- rowMeans(x.norm)
  message("calculating predictions...")
  t1 <- Sys.time()
  cvt <- expr.predict(t(x.est[-pred.genes[1], ]), x[pred.genes[1], ]/sf,
                      dfmax, nfolds)
  t2 <- Sys.time()
  t.diff <- t2-t1
  units(t.diff) <- "secs"
  if (parallel) {
    nworkers <- foreach::getDoParWorkers()
    t3 <- t.diff*1.1*npred/nworkers + 0.01*ngenes
    units(t3) <- "mins"
    message("Approximate finish time: ", t2+t3)
    if (nworkers == 1) {
      message("Only one worker assigned! Running sequentially...")
      mu.par <- matrix(0, npred, ncells)
      for (i in 1:length(pred.genes)) {
        mu.par[i, ] <- expr.predict(t(x.est[-pred.genes[i], ]),
                                    x[pred.genes[i], ]/sf, dfmax, nfolds,
                                    seed = pred.genes[i])
      }
    } else {
      message("Running in parallel: ", nworkers, " workers")
      gene.list <- chunk2(pred.genes, nworkers)
      mu.par <- foreach::foreach(i = 1:nworkers, .combine = rbind,
                                 .packages = c("glmnet", "SAVER")) %dopar% {
        mu.temp <- matrix(0, length(gene.list[[i]]), ncells)
        for (j in 1:length(gene.list[[i]])) {
          mu.temp[j, ] <- expr.predict(t(x.est[-gene.list[[i]][j], ]),
                                       x[gene.list[[i]][j], ]/sf, dfmax,
                                       nfolds, seed = gene.list[[i]][j])
        }
        return(mu.temp)
      }
    }

  } else {
    t3 <- t.diff*npred + 0.01*ngenes
    units(t3) <- "mins"
    message("Approximate finish time: ", t2+t3)
    mu.par <- matrix(0, npred, ncells)
    for (i in 1:length(pred.genes)) {
      mu.par[i, ] <- expr.predict(t(x.est[-pred.genes[i], ]),
                                  x[pred.genes[i], ]/sf, dfmax, nfolds,
                                  seed = pred.genes[i])
    }
  }
  message("Predictions finished. Calculating posterior...")
  if (pred.genes.only) {
    out <- lapply(1:4, function(x) matrix(0, npred, ncells))
    mu <- mu.par
    for (i in 1:length(pred.genes)) {
      post <- calc.post(x[pred.genes[i], ], mu[i, ], sf, scale.sf)
      out[[1]][i, ] <- post$estimate
      out[[2]][i, ] <- post$alpha
      out[[3]][i, ] <- post$beta
    }
    out[[4]] <- mu
    gene.names <- rownames(x)[pred.genes]
    cell.names <- colnames(x)
  } else {
    out <- lapply(1:4, function(x) matrix(0, ngenes, ncells))
    mu <- matrix(gene.means, ngenes, ncells)
    mu[pred.genes, ] <- mu.par
    for (i in 1:ngenes) {
      post <- calc.post(x[i, ], mu[i, ], sf, scale.sf)
      out[[1]][i, ] <- post$estimate
      out[[2]][i, ] <- post$alpha
      out[[3]][i, ] <- post$beta
    }
    out[[4]] <- mu
    gene.names <- rownames(x)
    cell.names <- colnames(x)
  }
  out.named <- lapply(out, function(x) {rownames(x) <- gene.names;
                                        colnames(x) <- cell.names; x})
  names(out.named) <- c("estimate", "alpha", "beta", "predicted")
  message("Done!")
  return(out.named)
}


