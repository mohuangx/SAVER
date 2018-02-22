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
#' @param do.fast Approximates the prediction step. Default is TRUE.
#'
#' @param size.factor Vector of cell size normalization factors.
#' If \code{x} is already normalized or normalization is not desired, use
#' \code{size.factor = 1}. Default uses mean library size normalization.
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
#' @param verbose If TRUE, prints index of gene
#'
#' @param predict.time If TRUE, calculates approximate finish time.
#'
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression}
#' \item{\code{se}}{Standard error of estimates}
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
#' system.time(saver(linnarsson, npred = 20))
#' registerDoParallel(cores = 2)
#' system.time(saver(linnarsson, npred = 20))
#' }
#'
#' @export
saver <- function(x, do.fast = TRUE, size.factor = NULL, npred = NULL, 
                  pred.cells = NULL, pred.genes = NULL, 
                  pred.genes.only = FALSE, null.model = FALSE, 
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
  pred.genes <- get.pred.genes(x, pred.genes, npred, ngenes)
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

  out <- saver.fit(x, x.est, do.fast, sf, scale.sf, pred.genes, pred.cells, 
                   pred.genes.only, null.model, verbose, predict.time, 
                   ngenes = nrow(x), ncells = ncol(x), gene.names, 
                   cell.names)
  class(out) <- "saver"
  message("Done!")
  out
}


