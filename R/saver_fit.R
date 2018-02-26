#' #' Fits SAVER
#' #'
#' #' Fits SAVER object
#' #'
#' #' The SAVER method starts by estimating the prior mean and variance for the
#' #' true expression level for each gene and cell. The prior mean is obtained
#' #' through predictions from a LASSO Poisson regression for each gene
#' #' implemented using the \code{glmnet} package. Then, the variance is estimated
#' #' through maximum likelihood assuming constant variance, Fano factor, or
#' #' coefficient of variation variance structure for each gene. The posterior
#' #' distribution is calculated and the posterior mean is reported as the SAVER
#' #' estimate.
#' #'
#' #' @param x An expression count matrix. The rows correspond to genes and
#' #' the columns correspond to cells.
#' #'
#' #' @param x.est The log-normalized predictor matrix. The rows correspond to
#' #' cells and the columns correspond to genes.
#' #'
#' #' @param do.fast Approximates the prediction step. Default is TRUE.
#' #'
#' #' @param sf Normalized size factor.
#' #'
#' #' @param scale.sf Scale of size factor.
#' #'
#' #' @param pred.genes Index of genes to perform regression prediction.
#' #'
#' #' @param pred.cells Index of cells to perform regression prediction.
#' #'
#' #' @param null.model Whether to use mean gene expression as prediction.
#' #'
#' #' @param ngenes Number of genes.
#' #'
#' #' @param ncells Number of cells.
#' #'
#' #' @param gene.names Name of genes.
#' #'
#' #' @param cell.names Name of cells.
#' #'
#' #' @return A list with the following components
#' #' \item{\code{estimate}}{Recovered (normalized) expression}
#' #' \item{\code{se}}{Standard error of estimates}
#' #' \item{\code{info}}{Information about fit}
#' #'
#' #' @export
#' #'
#' 
#' saver.fit <- function(x, x.est, do.fast, sf, scale.sf, pred.genes, pred.cells,
#'                       null.model, ngenes = nrow(x), ncells = ncol(x),
#'                       gene.names = rownames(x), cell.names = colnames(x)) {
#'   est <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
#'   se <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
#'   info <- c(list(0), rep(list(rep(0, ngenes)), 6), list(0), list(0))
#'   names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min",
#'                    "sd.cv", "pred.time", "var.time", "cutoff", "total.time")
#'   info$size.factor <- scale.sf*sf
#' 
#'   nworkers <- foreach::getDoParWorkers()
#'   message("Running SAVER with ", nworkers, " worker(s)")
#' 
#'   npred <- length(pred.genes)
#' 
#'   if (!null.model) {
#'     message("Calculating predictions for ", npred,
#'             " genes using ", ncol(x.est), " genes and ", nrow(x.est),
#'             " cells...")
#'   } else {
#'     message("Using means as predictions.")
#'   }
#'   set.seed(1)
#'   st <- proc.time()
#'   if (npred < ngenes) {
#'     ind <- c(sample(pred.genes, npred), sample((1:ngenes)[-pred.genes],
#'                                                ngenes-npred))
#'   } else {
#'     ind <- sample(1:ngenes, ngenes)
#'   }
#'   if (do.fast & !null.model) {
#'     if (npred < 100) {
#'       out1 <- calc.estimate(x[ind, ], x.est, cutoff = 0, coefs = NULL, sf,
#'                             scale.sf, gene.names[pred.genes], pred.cells,
#'                             null.model, nworkers, calc.maxcor = FALSE)
#'       est[ind, ] <- out1$est
#'       se[ind, ] <- out1$se
#'       for (j in 1:6) {
#'         info[[j+1]][ind] <- out1[[j+2]]
#'       }
#'     } else {
#'       ind1 <- ind[1:min(100, ngenes)]
#'       ptc <- proc.time()
#'       message("Step 1:")
#' 
#'       out1 <- calc.estimate(x[ind1, ], x.est, cutoff = 0, coefs = NULL, sf,
#'                             scale.sf, gene.names[pred.genes], pred.cells,
#'                             null.model, nworkers, calc.maxcor = TRUE)
#'       fit <- lm(sqrt(out1$sd.cv) ~ out1$maxcor)
#' 
#'       cutoff <- (0.5 - fit$coefficients[1])/fit$coefficients[2]
#'       info$cutoff <- unname(cutoff)
#' 
#'       est[ind1, ] <- out1$est
#'       se[ind1, ] <- out1$se
#'       for (j in 1:6) {
#'         info[[j+1]][ind1] <- out1[[j+2]]
#'       }
#' 
#'       n2 <- min(ceiling(200/mean(out1$maxcor > cutoff))+1, ngenes)
#' 
#'       if (npred < n2) {
#'         message("Step 2:")
#'         ind2 <- ind[101:npred]
#'         out2 <- calc.estimate(x[ind2, ], x.est, cutoff = 0, coefs = NULL, sf,
#'                               scale.sf, gene.names[pred.genes], pred.cells,
#'                               null.model, nworkers, calc.maxcor = TRUE)
#'         est[ind2, ] <- out2$est
#'         se[ind2, ] <- out2$se
#'         for (j in 1:6) {
#'           info[[j+1]][ind2] <- out2[[j+2]]
#'         }
#'         message("Step 3:")
#'         ind3 <- ind[(npred+1):ngenes]
#'         out3 <- calc.estimate(x[ind3, ], x.est, cutoff = 0, coefs = NULL, sf,
#'                               scale.sf, gene.names[pred.genes], pred.cells,
#'                               null.model, nworkers, calc.maxcor = TRUE)
#'         est[ind3, ] <- out3$est
#'         se[ind3, ] <- out3$se
#'         for (j in 1:6) {
#'           info[[j+1]][ind3] <- out3[[j+2]]
#'         }
#' 
#'       } else {
#'         message("Step 2:")
#'         ind2 <- ind[101:min(n2, length(ind))]
#'         out2 <- calc.estimate(x[ind2, ], x.est, cutoff = 0, coefs = NULL, sf,
#'                               scale.sf, gene.names[pred.genes], pred.cells,
#'                               null.model, nworkers, calc.maxcor = TRUE)
#'         est[ind2, ] <- out2$est
#'         se[ind2, ] <- out2$se
#'         for (j in 1:6) {
#'           info[[j+1]][ind2] <- out2[[j+2]]
#'         }
#'         maxcor <- c(out1$maxcor, out2$maxcor)
#'         pred <- which(maxcor > cutoff)
#'         lambda.max <- c(out1$lambda.max, out2$lambda.max)[pred]
#'         lambda.min <- c(out1$lambda.min, out2$lambda.min)[pred]
#'         coefs <- lm(log(lambda.max/lambda.min)^2 ~ maxcor[pred])$coefficients
#' 
#'         if (npred > n2) {
#'           ind3 <- ind[(n2+1):npred]
#' 
#'           ptc <- proc.time()
#'           message("Step 3:")
#' 
#'           out3 <- calc.estimate(x[ind3, ], x.est, cutoff, coefs, sf, scale.sf,
#'                                 gene.names[pred.genes], pred.cells, null.model,
#'                                 nworkers, calc.maxcor = TRUE)
#' 
#'           est[ind3, ] <- out3$est
#'           se[ind3, ] <- out3$se
#'           for (j in 1:6) {
#'             info[[j+1]][ind3] <- out3[[j+2]]
#'           }
#'         }
#'         if (npred < ngenes) {
#'           message("Step 4:")
#'           ind4 <- ind[(npred+1):length(ind)]
#'           out4 <- calc.estimate(x[ind4, ], x.est, cutoff, coefs, sf, scale.sf,
#'                                 gene.names[pred.genes], pred.cells,
#'                                 null.model = TRUE, nworkers,
#'                                 calc.maxcor = FALSE)
#' 
#'           est[ind4, ] <- out4$est
#'           se[ind4, ] <- out4$se
#'           for (j in 1:6) {
#'             info[[j+1]][ind4] <- out4[[j+2]]
#'           }
#'         }
#'       }
#'     }
#'   } else {
#'     out1 <- calc.estimate(x[ind, ], x.est, cutoff = 0, coefs = NULL, sf,
#'                           scale.sf, gene.names[pred.genes], pred.cells,
#'                           null.model, nworkers, calc.maxcor = FALSE)
#'     est[ind, ] <- out1$est
#'     se[ind, ] <- out1$se
#'     for (j in 1:6) {
#'       info[[j+1]][ind] <- out1[[j+2]]
#'     }
#'   }
#'   info[[8]] <- (proc.time - st)[3]
#'   list(estimate = est, se = se, info = info)
#' }
