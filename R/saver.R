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
#' @param parallel If \code{TRUE}, uses the parallel backend \code{foreach} to
#' parallelize prediction. Must register parallel backend beforehand
#' (\code{doParallel}, \code{doMC}). See example below.
#'
#' @param nzero Run prediction for genes that have at least this many
#' nonzero cells. Default is 10.
#'
#' @param npred Number of genes for regression prediction. Selects the top
#' \code{npred} genes in terms of mean expression for regression prediction.
#' Default is all genes.
#'
#' @param pred.cells Indices of cells to perform regression prediction.
#' Default is all cells.
#'
#' @param pred.genes Indices of specific genes to perform regression
#' prediction. Overrides \code{npred}. Default is all genes.
#'
#' @param pred.genes.only Return expression levels of only \code{pred.genes}.
#' Default is FALSE (returns expression levels of all genes).
#'
#' @param null.model Whether to use mean gene expression as prediction.
#'
#' @param dfmax Maximum number of genes used in prediction model. Default is
#' 300.
#'
#' @param nfolds Number of folds to be used in cross validation. Default is 5.
#'
#' @param remove.zero.genes Whether to remove genes with all zero.
#' Default is FALSE.
#'
#' @param verbose If TRUE, prints index of gene
#'
#' @param predict.time If TRUE, calculates approximate finish time.
#'
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression}
#' \item{\code{alpha}}{Posterior Gamma shape parameter}
#' \item{\code{beta}}{Posterior Gamma rate parameter}
#' \item{\code{info}}{Information about dataset}
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
saver <- function(x, size.factor = NULL, parallel = FALSE, nzero = 10,
                  npred = NULL, pred.cells = NULL, pred.genes = NULL,
                  pred.genes.only = FALSE, null.model = FALSE, dfmax = 300,
                  nfolds = 5, remove.zero.genes = FALSE, verbose = FALSE,
                  predict.time = TRUE) {
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  if (min(colSums(x)) == 0) {
    nzerocells <- sum(colSums(x) == 0)
    x <- x[, colSums(x) != 0]
    message("Removing ", nzerocells, " cells with zero expression.")
  }
  nonzero <- which(rowSums(x) != 0)
  if (remove.zero.genes) {
    x <- x[nonzero, ]
    message("Removing ", length(nonzero), " genes with zero expression.")
  }
  np <- dim(x)
  ngenes <- as.integer(np[1])
  ncells <- as.integer(np[2])
  if (verbose) {
    message(ngenes, " genes")
    message(ncells, " cells")
  }
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
  if (!is.null(pred.cells)) {
    if (!is.integer(pred.cells) | min(pred.cells) < 1 |
        max(pred.cells) > ncells) {
      stop("pred.cells must be column indices of x")
    }
  } else {
    pred.cells <- 1:ncells
  }
  if (!is.null(pred.genes)) {
    if (!is.integer(pred.genes) | min(pred.genes) < 1 |
        max(pred.genes) > ngenes) {
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
  good.genes <- which(rowSums(x > 0) >= nzero)
  x.est <- t(log(sweep(x[good.genes, ] + 1, 2, sf, "/")))
  if (pred.genes.only) {
    genes <- pred.genes
  } else {
    genes <- 1:ngenes
  }
  lasso.genes <- intersect(good.genes, pred.genes)
  nonlasso.genes <- genes[!(genes %in% lasso.genes)]
  nvar.vec <- rep(0, ngenes)
  if (!null.model) {
    message("Calculating predictions for ", length(lasso.genes),
            " genes using ", length(pred.cells), " cells...")
    if (npred > 5 & predict.time) {
      mu <- matrix(0, 5, ncells)
      s <- sample(1:length(good.genes), 5)
      t1 <- Sys.time()
      for (i in 1:5) {
        cvt <- system.time(mu[i, ] <- expr.predict(x.est[, -s[i]],
                                                   x[good.genes[s[i]], ]/sf,
                                                   pred.cells, dfmax,
                                                   nfolds)$mu)
        if (verbose) {
          print(cvt[3])
        }
      }
      t2 <- Sys.time()
      for (i in 1:5) {
        post <- calc.post(x[good.genes[s[i]], ], mu[i, ], sf, scale.sf)
      }
      t3 <- Sys.time()
      t.diff1 <- (t2-t1)/5
      t.diff2 <- (t3-t2)/5
      units(t.diff1) <- "secs"
      units(t.diff2) <- "secs"
      if (verbose) {
        message("Prediction: ", t.diff1, " per gene")
        message("Posterior calc: ", t.diff2, " per gene")
      }
    }
    nworkers <- foreach::getDoParWorkers()
    if (parallel & nworkers > 1) {
      if (npred > 5 & predict.time) {
        npred2 <- length(lasso.genes)
        npred3 <- length(nonlasso.genes)
        t4 <- t.diff1*npred2/nworkers*1.1 +
          t.diff2*(npred2/nworkers+npred3)*1.1
        units(t4) <- "mins"
        message("Approximate finish time: ", t3+t4)
      }
      message("Running in parallel: ", nworkers, " workers")
      lasso <- foreach::foreach(i = lasso.genes, .combine = "combine.mat",
                                .packages = c("glmnet", "SAVER")) %dopar% {
        ind <- which(i == good.genes)
        cv <- expr.predict(x.est[, -ind], x[i, ]/sf, pred.cells, dfmax, nfolds,
                           seed = i)
        mu <- cv$mu
        post <- calc.post(x[i, ], mu, sf, scale.sf)
        est <- post$estimate
        alpha <- post$alpha
        beta <- post$beta
        if (verbose) {
          message(i)
        }
        return(list(estimate = est, alpha = alpha, beta = beta,
                    nvar = cv$nvar))
      }
      out <- lapply(1:3, function(x) matrix(0, ngenes, ncells))
      for (i in 1:3) {
        out[[i]][lasso.genes, ] <- lasso[[i]]
      }
      nvar.vec[lasso.genes] <- lasso[[4]]
      if (length(nonlasso.genes) > 0) {
        nvar.vec[nonlasso.genes] <- 0
        for (i in nonlasso.genes) {
          post <- calc.post(x[i, ], mean(x[i, ]/sf), sf, scale.sf)
          out[[1]][i, ] <- post$estimate
          out[[2]][i, ] <- post$alpha
          out[[3]][i, ] <- post$beta
        }
      }
    } else {
      if (parallel & nworkers == 1) {
        message("Only one worker assigned! Running sequentially...")
      }
      if (npred > 5 & predict.time) {
        npred2 <- length(lasso.genes)
        t4 <- t.diff1*npred2*1.1 + t.diff2*length(genes)*1.1
        units(t4) <- "mins"
        message("Approximate finish time: ", t3+t4)
      }
      out <- lapply(1:3, function(x) matrix(0, ngenes, ncells))
      k <- 0
      for (j in genes) {
        k <- k+1
        if (j %in% lasso.genes) {
          ind <- which(j == good.genes)
          cv <- expr.predict(x.est[, -ind], x[j, ]/sf, pred.cells, dfmax,
                             nfolds, seed = j)
          mu <- cv$mu
          nvar <- cv$nvar
          if (verbose) {
            print(j)
          }
        } else {
          mu <- rep(mean(x[j, ]/sf), ncells)
          nvar <- 0
        }
        post <- calc.post(x[j, ], mu, sf, scale.sf)
        out[[1]][j, ] <- post$estimate
        out[[2]][j, ] <- post$alpha
        out[[3]][j, ] <- post$beta
        nvar.vec[j] <- nvar
      }
    }
  } else {
    message("Using means as predictions.")
    if (npred > 5 & predict.time) {
      mu <- matrix(0, 5, ncells)
      s <- sample(1:length(good.genes), 5)
      t2 <- Sys.time()
      for (i in 1:5) {
        mu[i, ] <- rep(mean(x[good.genes[s[i]], ]/sf), ncells)
        post <- calc.post(x[good.genes[s[i]], ], mu[i, ], sf, scale.sf)
      }
      t3 <- Sys.time()
      t.diff2 <- (t3-t2)/5
      units(t.diff2) <- "secs"
      if (verbose) {
        message("Posterior calc: ", t.diff2, " per gene")
      }
    }
    if (npred > 5 & predict.time) {
      npred2 <- length(lasso.genes)
      t4 <- t.diff2*length(genes)*1.1
      units(t4) <- "mins"
      message("Approximate finish time: ", t3+t4)
    }
    out <- lapply(1:3, function(x) matrix(0, ngenes, ncells))
    for (j in genes) {
      mu <- rep(mean(x[j, ]/sf), ncells)
      nvar <- 0
      post <- calc.post(x[j, ], mu, sf, scale.sf)
      out[[1]][j, ] <- post$estimate
      out[[2]][j, ] <- post$alpha
      out[[3]][j, ] <- post$beta
      nvar.vec[j] <- nvar
    }
  }
  if (pred.genes.only) {
    for (i in 1:3) {
      out[[i]] <- out[[i]][pred.genes, ]
    }
    nvar.vec <- nvar.vec[pred.genes]
    gene.names <- rownames(x)[pred.genes]
    cell.names <- colnames(x)
  } else {
    gene.names <- rownames(x)
    cell.names <- colnames(x)
  }
  out.named <- lapply(out[1:3], function(x) {rownames(x) <- gene.names;
  colnames(x) <- cell.names; x})
  info <- list(size.factor = sf*scale.sf, nvar = nvar.vec, ngenes = ngenes,
               ncells = ncells, genes = gene.names, cells = cell.names)
  out.named[[4]] <- info
  names(out.named[[4]][[1]]) <- cell.names
  names(out.named[[4]][[2]]) <- gene.names
  names(out.named) <- c("estimate", "alpha", "beta", "info")
  class(out.named) <- "saver"
  message("Done!")
  return(out.named)
}


