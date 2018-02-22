#' Fits SAVER
#'
#' Fits SAVER object
#'
#' Fits SAVER object
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
#' @export
#' 

saver.fit <- function(x, x.est, do.fast, sf, scale.sf, pred.genes, pred.cells, 
                      pred.genes.only, null.model, verbose, predict.time, 
                      ngenes = nrow(x), ncells = ncol(x),
                      gene.names = rownames(x), cell.names = colnames(x)) {
  est <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  se <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  info <- c(list(0), rep(list(rep(0, ngenes)), 6), list(0))
  names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min", 
                   "sd.cv", "pred.time", "var.time", "cutoff")
  info$size.factor <- scale.sf*sf
  
  nworkers <- foreach::getDoParWorkers()
  message("Running SAVER with ", nworkers, " worker(s)")
  
  npred <- length(pred.genes)
  
  if (!null.model) {
    message("Calculating predictions for ", ngenes,
            " genes using ", ncol(x.est), " genes and ", nrow(x.est),
            " cells...")
  } else {
    message("Using means as predictions.")
  }
  set.seed(1)
  ind <- sample(1:ngenes, ngenes)
  if (do.fast & !null.model) {
    if (npred < 100) {
      out1 <- calc.estimate(x[ind, ], x.est, cutoff = 0, coefs = NULL, sf, 
                            scale.sf, gene.names[pred.genes], pred.cells, 
                            null.model, nworkers, verbose)
      est <- out1$est
      se <- out1$est
      for (j in 1:6) {
        info[[j+1]] <- out1[[j+2]]
      }
    } else {
      ind1 <- ind[1:min(100, length(ind))]
      ptc <- proc.time()
      message("Calculating cutoff for ", length(ind1), " genes")
      
      out1 <- calc.estimate(x[ind1, ], x.est, cutoff = 0, coefs = NULL, sf, 
                            scale.sf, gene.names[pred.genes], pred.cells, 
                            null.model, nworkers, verbose)
      fit <- lm(sqrt(out1$sd.cv) ~ out1$maxcor)
      
      cutoff <- (0.5 - fit$coefficients[1])/fit$coefficients[2]
      info$cutoff <- unname(cutoff)
      
      est[ind1, ] <- out1$est
      se[ind1, ] <- out1$se
      for (j in 1:6) {
        info[[j+1]][ind1] <- out1[[j+2]]
      }
      message((proc.time()-ptc)[3])
      message("cutoff is: ", cutoff)
      message("proportion of genes over cutoff: ", mean(out1$maxcor > cutoff))
      
      if (npred < n2)
      
      n2 <- ceiling(200/mean(out1$maxcor > cutoff))+1
      ind2 <- ind[101:min(n2, length(ind))]
    }
    
    
    if (length(genes) > 100) {
      n2 <- ceiling(200/mean(out1$maxcor > cutoff))+1
      ind2 <- ind[101:min(n2, length(ind))]
      
      ptc <- proc.time()
      message("cutoff is: ", cutoff)
      message("proportion of genes over cutoff: ", mean(out1$maxcor > cutoff))
      message("Calculating fit for ", length(ind2), " genes")
      
      out2 <- calc.estimate(x[ind2, ], x.est, cutoff, coefs = NULL, sf, 
                            scale.sf, gene.names[pred.genes], pred.cells, 
                            null.model, nworkers, verbose)
      
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
      coefs <- lm(log(lambda.max/lambda.min)^2 ~ maxcor[pred])$coefficients
      
      ind3 <- ind[(n2+1):length(ind)]
      
      ptc <- proc.time()
      message("Calculating rest for ", length(ind3), " genes")
      
      out3 <- calc.estimate(x[ind3, ], x.est, cutoff, coefs, sf, scale.sf, 
                            gene.names[pred.genes], pred.cells, null.model,
                            nworkers, verbose)
      
      est[ind3, ] <- out3$est
      se[ind3, ] <- out3$se
      for (j in 1:6) {
        info[[j+1]][ind3] <- out3[[j+2]]
      }
    }
    message((proc.time()-ptc)[3])
    
  } else {
    out1 <- calc.estimate(x[ind, ], x.est, cutoff = 0, coefs = NULL, sf, 
                          scale.sf, gene.names[pred.genes], pred.cells, 
                          null.model, nworkers, verbose)
    est <- out1$est
    se <- out1$est
    for (j in 1:6) {
      info[[j+1]] <- out1[[j+2]]
    }
  }
  list(estimate = est, se = se, info = info)
}
