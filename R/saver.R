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
#' the columns correspond to cells. Can be sparse.
#'
#' @param do.fast Approximates the prediction step. Default is TRUE.
#'
#' @param ncores Number of cores to use. Default is 1.
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
#' @param mu Matrix of prior means.
#'
#' @param estimates.only Only return SAVER estimates. Default is FALSE.
#'
#'
#' @return If `estimates.only = TRUE`, then a matrix of SAVER estimates.
#'
#' If `estimates.only = FALSE`, a list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression.}
#' \item{\code{se}}{Standard error of estimates.}
#' \item{\code{info}}{Information about dataset.}
#'
#' The \code{info} element is a list with the following components:
#' \item{\code{size.factor}}{Size factor used for normalization.}
#' \item{\code{maxcor}}{Maximum absolute correlation for each gene. 2 if not
#' calculated}
#' \item{\code{lambda.max}}{Smallest value of lambda which gives the null
#' model.}
#' \item{\code{lambda.min}}{Value of lambda from which the prediction model is
#' used}
#' \item{\code{sd.cv}}{Difference in the number of standard deviations in
#' deviance between the model with lowest cross-validation error and the null
#' model}
#' \item{\code{pred.time}}{Time taken to generate predictions.}
#' \item{\code{var.time}}{Time taken to estimate variance.}
#' \item{\code{maxcor}}{Maximum absolute correlation cutoff used to determine
#' if a gene should be predicted.}
#' \item{\code{lambda.coefs}}{Coefficients for estimating lambda with lowest
#' cross-validation error.}
#' \item{\code{total.time}}{Total time for SAVER estimation.}
#'
#' @examples
#' data("linnarsson")
#' 
#' \dontrun{
#' system.time(linnarsson_saver <- saver(linnarsson, ncores = 12))
#' }
#'
#' # predictions for top 5 highly expressed genes
#' \dontrun{
#' saver2 <- saver(linnarsson, npred = 5)
#' }
#'
#' # predictions for certain genes
#' \dontrun{
#' genes <- c("Thy1", "Mbp", "Stim2", "Psmc6", "Rps19")
#' genes.ind <- which(rownames(linnarsson) %in% genes)
#' saver3 <- saver(linnarsson, pred.genes = genes.ind)
#' }
#'
#' # return only certain genes
#' \dontrun{
#' saver4 <- saver(linnarsson, pred.genes = genes.ind, pred.genes.only = TRUE)
#' }
#'
#' @import foreach
#'
#' @export
saver <- function(x, do.fast = TRUE, ncores = 1, size.factor = NULL,
                  npred = NULL, pred.cells = NULL, pred.genes = NULL,
                  pred.genes.only = FALSE, null.model = FALSE, mu = NULL,
                  estimates.only = FALSE) {
  if (!is.null(mu)) {
    mu <- check.mu(x, mu)
  }
  x <- clean.data(x)
  np <- dim(x)
  ngenes <- as.integer(np[1])
  ncells <- as.integer(np[2])

  message(ngenes, " genes, ", ncells, " cells")

  # assign size factor
  sf.out <- calc.size.factor(x, size.factor, ncells)
  sf <- sf.out[[1]]
  scale.sf <- sf.out[[2]]

  gene.names <- rownames(x)
  cell.names <- colnames(x)

  # Set up parallel backend
  cl.create <- FALSE
  if (foreach::getDoParWorkers() > 1) {
    if (ncores > 1 & foreach::getDoParWorkers() != ncores) {
      message(paste("Parallel backend already registered and is inconsistent",
                    "with ncores."))
    }
    ncores <- getDoParWorkers()
  } else {
    if (ncores > 1) {
      cl.create <- TRUE
      cl <- parallel::makeCluster(ncores, outfile = "")
      doParallel::registerDoParallel(cl)
      on.exit({
        parallel::stopCluster(cl)
        foreach::registerDoSEQ()
      })
    }
  }

  # if prior means are provided
  if (!is.null(mu)) {
    out <- saver.fit.mean(x, ncores, sf, scale.sf, mu, ngenes = nrow(x),
                          ncells = ncol(x), gene.names, cell.names,
                          estimates.only)
  } else if (null.model) {
    out <- saver.fit.null(x, ncores, sf, scale.sf, ngenes = nrow(x),
                          ncells = ncol(x), gene.names, cell.names,
                          estimates.only)
  } else {
    # assign pred.cells and pred.genes
    pred.cells <- get.pred.cells(pred.cells, ncells)
    pred.genes <- get.pred.genes(x, pred.genes, npred, ngenes)
    npred <- length(pred.genes)


    good.genes <- which(Matrix::rowMeans(sweep(x, 2, sf, "/")) >= 0.1)
    x.est <- t(as.matrix(log(sweep(x[good.genes, ] + 1, 2, sf, "/"))))
    if (pred.genes.only) {
      x <- x[pred.genes, , drop = FALSE]
      pred.genes <- 1:nrow(x)
      gene.names <- rownames(x)
      cell.names <- colnames(x)
    }

    out <- saver.fit(x, x.est, do.fast, ncores, sf, scale.sf, pred.genes,
                     pred.cells, null.model, ngenes = nrow(x),
                     ncells = ncol(x), gene.names, cell.names,
                     estimates.only)
  }
  message("Done!")
  message("Finish time: ", Sys.time())
  message("Total time: ", format(out$info$total.time))
  if (!estimates.only) {
    class(out) <- "saver"
    out
  } else {
    out$estimate
  }
}


