#' Runs SAVER-fast
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
#' @param nlambda Number of lambda to calculate in cross validation. Default is
#' 5.
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
saver_fast <- function(x, size.factor = NULL, parallel = FALSE, nzero = 10,
                  npred = NULL, pred.cells = NULL, pred.genes = NULL,
                  output.se = TRUE, 
                  pred.genes.only = FALSE, null.model = FALSE, dfmax = 300,
                  nfolds = 5, nlambda = 50, remove.zero.genes = FALSE,
                  verbose = FALSE, predict.time = TRUE) {
  x <- clean.data(x)
  np <- dim(x)
  ngenes <- as.integer(np[1])
  ncells <- as.integer(np[2])
  if (verbose) {
    message(ngenes, " genes")
    message(ncells, " cells")
  }
  
  # assign size factor
  sf.out <- calc.size.factor(x, size.factor, ncells)
  sf <- sf.out[[1]]
  scale.sf <- sf.out[[2]]
  
  # assign pred.cells and pred.genes
  pred.cells <- get.pred.cells(pred.cells, ncells)
  pred.genes <- get.pred.genes(pred.genes, npred, ngenes)
  npred <- length(pred.genes)
  
  good.genes <- intersect(which(rowMeans(sweep(x, 2, sf, "/")) >= 0.1),
                          pred.genes)
  x.est <- t(log(sweep(x[good.genes, ] + 1, 2, sf, "/")))
  if (pred.genes.only) {
    genes <- pred.genes
  } else {
    genes <- 1:ngenes
  }
  x <- x[genes, ]
  est <- matrix(0, length(genes), ncells)
  se <- matrix(0, length(genes), ncells)
  
  nworkers <- foreach::getDoParWorkers()
  message("Running SAVER with ", nworkers, " worker(s)")
  
  if (!null.model) {
    message("Calculating predictions for ", ngenes,
            " genes using ", ncol(x.est), " genes and ", nrow(x.est),
            " cells...")
  } else {
    message("Using means as predictions.")
  }
  set.seed(1)
  ind <- sample(genes, length(genes))
  ind1 <- ind[1:min(100, length(ind))]
  
  out1 <- calc.cutoff(x[ind1, ], x.est, sf, npred, pred.cells, nworkers, 
                      output.se, verbose, index = NULL)
  
  est[ind1, ] <- out1$est
  se[ind1, ] <- out1$se
  
  if (length(genes) > 100) {
    cutoff <- out1$cutoff
    n2 <- ceiling(200/mean(out1$maxcor > cutoff))+1
    ind2 <- ind[101:min(n2, length(ind))]
    out2 <- calc.lambda(x[ind2, ], x.est, cutoff, sf, npred, pred.cells, 
                        nworkers, output.se, verbose)
    
    est[ind2, ] <- out2$est
    se[ind2, ] <- out2$se
  }
  
  if (length(genes) > n2) {
    maxcor <- c(out1$maxcor, out2$maxcor)
    pred <- which(maxcor > cutoff)
    lambda.max <- c(out1$lambda.max, out2$lambda.max)[pred]
    lambda.min <- c(out1$lambda.min, out2$lambda.min)[pred]
    fit <- lm(log(lambda.max/lambda.min)^2 ~ maxcor[pred])
    
    ind3 <- ind[(n2+1):length(ind)]
    
    out3 <- calc.estimate(x[ind3[1:100], ], x.est, cutoff, fit, sf, npred, pred.cells,
                          nworkers, output.se, verbose)
    
    est[ind3, ] <- out3$est
    se[ind3, ] <- out3$se
  }
  
  

  lasso.genes <- intersect(good.genes, pred.genes)
  nonlasso.genes <- genes[!(genes %in% lasso.genes)]
  nvar.vec <- rep(0, ngenes)
  sd.vec <- rep(0, ngenes)
  if (!null.model) {

    if (npred > 5 & predict.time) {
      mu <- matrix(0, 5, ncells)
      s <- sample(1:length(good.genes), 5)
      t1 <- Sys.time()
      for (i in 1:5) {
        cvt <- system.time(mu[i, ] <- expr.predict(x.est[, -s[i]],
                                                   x[good.genes[s[i]], ]/sf,
                                                   pred.cells, dfmax,
                                                   nfolds, nlambda)$mu)
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
      lasso <- foreach::foreach(xit = iterators::iter(x[lasso.genes, ], by = "row"), 
                                i = lasso.genes,
                                .packages = c("glmnet", "SAVER",
                                              "bigmemory", "iterators")) %dopar% {
                                                ind <- which(i == good.genes)
                                                if (length(ind) == 0) {
                                                  cv <- expr.predict(x.est, c(xit)/sf, pred.cells, dfmax, nfolds,
                                                                     nlambda, seed = i)
                                                } else {
                                                  cv <- expr.predict(x.est[, -ind], c(xit)/sf, pred.cells, dfmax, nfolds,
                                                                     nlambda, seed = i)
                                                }
                                                mu <- cv$mu
                                                post <- calc.post(c(xit), mu, sf, scale.sf)
                                                est <- unname(post$estimate)
                                                alpha <- unname(post$alpha)
                                                beta <- unname(post$beta)
                                                if (verbose) {
                                                  message(i)
                                                }
                                                return(list(estimate = est, alpha = alpha, beta = beta,
                                                            nvar = unname(cv$nvar),
                                                            sd.cv = unname(cv$sd.cv)))
                                              }
      out <- lapply(1:3, function(x) matrix(0, ngenes, ncells))
      for (i in 1:3) {
        tempvec <- lapply(lasso, `[[`, i)
        out[[i]][lasso.genes, ] <- matrix(unlist(tempvec),
                                          nrow = length(tempvec), byrow = TRUE)
      }
      nvar.vec[lasso.genes] <- unlist(lapply(lasso, `[[`, 4))
      sd.vec[lasso.genes] <- unlist(lapply(lasso, `[[`, 5))
      if (length(nonlasso.genes) > 0) {
        nvar.vec[nonlasso.genes] <- 0
        sd.vec[nonlasso.genes] <- 0
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
          if (length(ind) == 0) {
            cv <- expr.predict(x.est, x[j, ]/sf, pred.cells, dfmax, nfolds,
                               nlambda, seed = j)
          } else {
            cv <- expr.predict(x.est[, -ind], x[j, ]/sf, pred.cells, dfmax, nfolds,
                               nlambda, seed = j)
          }
          mu <- cv$mu
          nvar <- cv$nvar
          sd.cv <- cv$sd.cv
          if (verbose) {
            print(j)
          }
        } else {
          mu <- rep(mean(x[j, ]/sf), ncells)
          nvar <- 0
          sd.cv <- 0
        }
        post <- calc.post(x[j, ], mu, sf, scale.sf)
        out[[1]][j, ] <- post$estimate
        out[[2]][j, ] <- post$alpha
        out[[3]][j, ] <- post$beta
        nvar.vec[j] <- nvar
        sd.vec[j] <- sd.cv
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
      sd.cv <- 0
      post <- calc.post(x[j, ], mu, sf, scale.sf)
      out[[1]][j, ] <- post$estimate
      out[[2]][j, ] <- post$alpha
      out[[3]][j, ] <- post$beta
      nvar.vec[j] <- nvar
      sd.vec[j] <- sd.cv
    }
  }
  if (pred.genes.only) {
    for (i in 1:3) {
      out[[i]] <- out[[i]][pred.genes, , drop = FALSE]
    }
    nvar.vec <- nvar.vec[pred.genes]
    sd.vec <- sd.vec[pred.genes]
    gene.names <- rownames(x)[pred.genes]
    cell.names <- colnames(x)
  } else {
    gene.names <- rownames(x)
    cell.names <- colnames(x)
  }
  out.named <- lapply(out[1:3], function(x) {rownames(x) <- gene.names;
  colnames(x) <- cell.names; x})
  info <- list(size.factor = sf*scale.sf, nvar = nvar.vec, sd.cv = sd.vec,
               ngenes = ngenes, ncells = ncells, genes = gene.names,
               cells = cell.names)
  out.named[[4]] <- info
  names(out.named[[4]][[1]]) <- cell.names
  names(out.named[[4]][[2]]) <- gene.names
  names(out.named[[4]][[3]]) <- gene.names
  names(out.named) <- c("estimate", "alpha", "beta", "info")
  class(out.named) <- "saver"
  message("Done!")
  return(out.named)
}


