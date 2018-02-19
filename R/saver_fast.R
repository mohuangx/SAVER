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
  
  gene.names <- rownames(x)
  cell.names <- colnames(x)
  
  good.genes <- intersect(which(rowMeans(sweep(x, 2, sf, "/")) >= 0.1),
                          pred.genes)
  x.est <- t(log(sweep(x[good.genes, ] + 1, 2, sf, "/")))
  if (pred.genes.only) {
    genes <- pred.genes
  } else {
    genes <- 1:ngenes
  }
  x <- x[genes, ]
  est <- matrix(0, length(genes), ncells, dimnames = list(gene.names[genes], 
                                                          cell.names))
  se <- matrix(0, length(genes), ncells, dimnames = list(gene.names[genes],
                                                         cell.names))
  info <- c(list(0), rep(list(rep(0, ngenes)), 6), list(0))
  names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min", 
                   "sd.cv", "pred.time", "var.time", "cutoff")
  info$size.factor <- scale.sf*sf
  
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
  ptc <- proc.time()
  message("Calculating cutoff for ", length(ind1), " genes")
  
  out1 <- calc.cutoff(x[ind1, ], x.est, sf, scale.sf, npred, pred.cells, 
                      nworkers, output.se, verbose, index = NULL)
  
  fit <- lm(sqrt(out1$sd.cv) ~ out1$maxcor)
  cutoff <- (0.5 - fit$coefficients[1])/fit$coefficients[2]
  info$cutoff <- unname(cutoff)
  
  est[ind1, ] <- out1$est
  se[ind1, ] <- out1$se
  for (j in 1:6) {
    info[[j+1]][ind1] <- out1[[j+2]]
  }
  message((proc.time()-ptc)[3])
  

  if (length(genes) > 100) {
    n2 <- ceiling(200/mean(out1$maxcor > cutoff))+1
    ind2 <- ind[101:min(n2, length(ind))]
    
    ptc <- proc.time()
    message("cutoff is: ", cutoff)
    message("proportion of genes over cutoff: ", mean(out1$maxcor > cutoff))
    message("Calculating fit for ", length(ind2), " genes")
    
    out2 <- calc.lambda(x[ind2, ], x.est, cutoff, sf, scale.sf, npred, 
                        pred.cells, nworkers, output.se, verbose)
    
    est[ind2, ] <- out2$est
    se[ind2, ] <- out2$se
    for (j in 1:6) {
      info[[j+1]][ind2] <- out2[[j+2]]
    }
    message((proc.time()-ptc)[3])
  }

  

  if (length(genes) > n2) {
    maxcor <- c(out1$maxcor, out2$maxcor)
    pred <- which(maxcor > cutoff)
    lambda.max <- c(out1$lambda.max, out2$lambda.max)[pred]
    lambda.min <- c(out1$lambda.min, out2$lambda.min)[pred]
    fit <- lm(log(lambda.max/lambda.min)^2 ~ maxcor[pred])
    
    ind3 <- ind[(n2+1):length(ind)]
    
    ptc <- proc.time()
    message("Calculating rest for ", length(ind3), " genes")
    
    out3 <- calc.estimate(x[ind3, ], x.est, cutoff, fit, sf, scale.sf, npred, 
                          pred.cells, nworkers, output.se, verbose)
    
    est[ind3, ] <- out3$est
    se[ind3, ] <- out3$se
    for (j in 1:6) {
      info[[j+1]][ind3] <- out3[[j+2]]
    }
  }
  message((proc.time()-ptc)[3])
  
  out <- list(estimate = est, se = se, info = info)
  class(out) <- "saver"
  message("Done!")
  out
}


