#' Fits SAVER
#'
#' Fits SAVER object
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
#' @param x.est The log-normalized predictor matrix. The rows correspond to
#' cells and the columns correspond to genes.
#'
#' @param do.fast Approximates the prediction step. Default is TRUE.
#'
#' @param sf Normalized size factor.
#'
#' @param scale.sf Scale of size factor.
#'
#' @param pred.genes Index of genes to perform regression prediction.
#'
#' @param pred.cells Index of cells to perform regression prediction.
#'
#' @param null.model Whether to use mean gene expression as prediction.
#'
#' @param ngenes Number of genes.
#'
#' @param ncells Number of cells.
#'
#' @param gene.names Name of genes.
#'
#' @param cell.names Name of cells.
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression}
#' \item{\code{se}}{Standard error of estimates}
#' \item{\code{info}}{Information about fit}
#'
#' @export
#'

saver.fit <- function(x, x.est, do.fast, sf, scale.sf, pred.genes, pred.cells,
                      null.model, ngenes = nrow(x), ncells = ncol(x),
                      gene.names = rownames(x), cell.names = colnames(x)) {
  est <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  se <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  info <- c(list(0), rep(list(rep(0, ngenes)), 6), list(0), list(0))
  names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min",
                   "sd.cv", "pred.time", "var.time", "cutoff", "total.time")
  info$size.factor <- scale.sf*sf
  
  nworkers <- foreach::getDoParWorkers()
  message("Running SAVER with ", nworkers, " worker(s)")
  
  npred <- length(pred.genes)
  message(npred)
  if (!null.model) {
    message("Calculating predictions for ", npred,
            " genes using ", ncol(x.est), " genes and ", nrow(x.est),
            " cells...")
  } else {
    message("Using means as predictions.")
  }
  set.seed(1)
  st <- Sys.time()
  if (npred < ngenes) {
    ind <- c(sample(pred.genes, npred), sample((1:ngenes)[-pred.genes],
                                               ngenes-npred))
  } else {
    ind <- sample(1:ngenes, ngenes)
  }
  if (do.fast & !null.model) {
    n1 <- min(max(8, nworkers), npred)
    ind1 <- ind[1:n1]
    message("Estimating finish time with ", length(ind1), " genes.")
    t1 <- Sys.time()
    out <- calc.estimate(x[ind1, , drop = FALSE], x.est, cutoff = 0,
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                         pred.cells, null.model, nworkers, calc.maxcor = TRUE)
    est[ind1, ] <- out$est
    se[ind1, ] <- out$se
    for (j in 1:6) {
      info[[j+1]][ind1] <- out[[j+2]]
    }
    
    if (n1 == npred) {
      if (n1 == ngenes) {
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      } else {
        ind5 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind5), " genes.")
        out <- calc.estimate(x[ind5, , drop = FALSE], x.est, cutoff = 0, 
                             coefs = NULL, sf, scale.sf, 
                             gene.names[pred.genes], pred.cells, 
                             null.model = TRUE, nworkers, calc.maxcor = FALSE)
        est[ind5, ] <- out$est
        se[ind5, ] <- out$se
        for (j in 1:6) {
          info[[j+1]][ind5] <- out[[j+2]]
        }
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      }
    }
    
    n2 <- min(max(100, nworkers), npred)
    ind2 <- ind[(n1+1):n2]
    
    t2 <- Sys.time()
    tdiff <- (t2-t1)/n1*n2/0.5 + max(0, (t2-t1)/n1*(npred-n2*3)/20) +
      (ngenes-npred)*mean(out[[8]])
    units(tdiff) <- "secs"
    message("Approximate finish time: ", Sys.time() + tdiff)

    message("Calculating max cor cutoff with ", length(ind2), " genes.")
    out <- calc.estimate(x[ind2, , drop = FALSE], x.est, cutoff = 0, 
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes], 
                         pred.cells, null.model, nworkers, calc.maxcor = TRUE)
    est[ind2, ] <- out$est
    se[ind2, ] <- out$se
    for (j in 1:6) {
      info[[j+1]][ind2] <- out[[j+2]]
    }
    
    if (n2 == npred) {
      if (n2 == ngenes) {
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      } else {
        ind5 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind5), " genes.")
        out <- calc.estimate(x[ind5, , drop = FALSE], x.est, cutoff = 0, 
                             coefs = NULL, sf, scale.sf, 
                             gene.names[pred.genes], pred.cells, 
                             null.model = TRUE, nworkers, calc.maxcor = FALSE)
        est[ind5, ] <- out$est
        se[ind5, ] <- out$se
        for (j in 1:6) {
          info[[j+1]][ind5] <- out[[j+2]]
        }
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      }
    }
    fit <- lm(sqrt(info$sd.cv[ind[1:n2]]) ~ info$maxcor[ind[1:n2]])
    cutoff <- (0.5 - fit$coefficients[1])/fit$coefficients[2]
    info$cutoff <- unname(cutoff)
    
    message(cutoff)

    n3 <- min(max(ceiling(n2/mean(info$maxcor[ind[1:n2]] > cutoff) + n2),
                  nworkers + n2), npred)
    ind3 <- ind[(n2+1):n3]
    
    message(n3)
    
    message("Calculating lambda coefficients with ", length(ind3), " genes.")
    out <- calc.estimate(x[ind3, , drop = FALSE], x.est, cutoff, coefs = NULL, 
                         sf, scale.sf, gene.names[pred.genes], pred.cells,
                         null.model, nworkers, calc.maxcor = TRUE)
    est[ind3, ] <- out$est
    se[ind3, ] <- out$se
    for (j in 1:6) {
      info[[j+1]][ind3] <- out[[j+2]]
    }
    
    if (n3 == npred) {
      if (n3 == ngenes) {
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      } else {
        ind5 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind5), " genes.")
        out <- calc.estimate(x[ind5, , drop = FALSE], x.est, cutoff = 0, 
                             coefs = NULL, sf, scale.sf, 
                             gene.names[pred.genes], pred.cells,
                             null.model = TRUE, nworkers, calc.maxcor = FALSE)
        est[ind5, ] <- out$est
        se[ind5, ] <- out$se
        for (j in 1:6) {
          info[[j+1]][ind5] <- out[[j+2]]
        }
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      }
    }
    
    pred <- which(info$maxcor > cutoff)
    lambda.max <- info$lambda.max[pred]
    lambda.min <- info$lambda.min[pred]
    coefs <- lm(log(lambda.max/lambda.min)^2 ~ info$maxcor[pred])$coefficients
    
    n4 <- npred
    ind4 <- ind[(n3+1):n4]
    
    message("Predicting ", length(ind4), " genes.")
    out <- calc.estimate(x[ind4, , drop = FALSE], x.est, cutoff, coefs, sf, 
                         scale.sf, gene.names[pred.genes], pred.cells, 
                         null.model, nworkers, calc.maxcor = TRUE)
    
    est[ind4, ] <- out$est
    se[ind4, ] <- out$se
    for (j in 1:6) {
      info[[j+1]][ind4] <- out[[j+2]]
    }
    
    if (n4 == ngenes) {
      info[[9]] <- as.numeric(Sys.time() - st)
      return(list(estimate = est, se = se, info = info))
    } else {
      ind5 <- ind[(npred+1):ngenes]
      message("Estimating remaining ", length(ind5), " genes.")
      out <- calc.estimate(x[ind5, , drop = FALSE], x.est, cutoff = 0, 
                           coefs = NULL, sf, scale.sf, gene.names[pred.genes], 
                           pred.cells, null.model = TRUE, nworkers, 
                           calc.maxcor = FALSE)
      est[ind5, ] <- out$est
      se[ind5, ] <- out$se
      for (j in 1:6) {
        info[[j+1]][ind5] <- out[[j+2]]
      }
      info[[9]] <- as.numeric(Sys.time() - st)
      return(list(estimate = est, se = se, info = info))
    }
  } else {
    n1 <- min(max(8, nworkers), npred)
    ind1 <- ind[1:n1]
    message("Estimating finish time with ", length(ind1), " genes.")
    t1 <- Sys.time()
    out <- calc.estimate(x[ind1, , drop = FALSE], x.est, cutoff = 0, 
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes], 
                         pred.cells, null.model, nworkers, calc.maxcor = FALSE)
    est[ind1, ] <- out$est
    se[ind1, ] <- out$se
    for (j in 1:6) {
      info[[j+1]][ind1] <- out[[j+2]]
    }
    
    if (n1 == npred) {
      if (n1 == ngenes) {
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      } else {
        ind5 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind5), " genes.")
        out <- calc.estimate(x[ind5, , drop = FALSE], x.est, cutoff = 0, 
                             coefs = NULL, sf, scale.sf, 
                             gene.names[pred.genes], pred.cells,
                             null.model = TRUE, nworkers, calc.maxcor = FALSE)
        est[ind5, ] <- out$est
        se[ind5, ] <- out$se
        for (j in 1:6) {
          info[[j+1]][ind5] <- out[[j+2]]
        }
        info[[9]] <- as.numeric(Sys.time() - st)
        return(list(estimate = est, se = se, info = info))
      }
    }
    
    n2 <- ngenes
    ind2 <- ind[(n1+1):n2]
    
    t2 <- Sys.time()
    tdiff <- (t2-t1)/n1*n2
    units(tdiff) <- "secs"
    message("Approximate finish time: ", Sys.time() + tdiff)
    
    message("Estimating remaining ", length(ind2), " genes.")
    out <- calc.estimate(x[ind2, , drop = FALSE], x.est, cutoff = 0, 
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes], 
                         pred.cells, null.model, nworkers, calc.maxcor = FALSE)
    est[ind, ] <- out$est
    se[ind, ] <- out$se
    for (j in 1:6) {
      info[[j+1]][ind] <- out[[j+2]]
    }
    info[[9]] <- as.numeric(Sys.time() - st)
    return(list(estimate = est, se = se, info = info))
  }
  info[[9]] <- as.numeric(Sys.time() - st)
  return(list(estimate = est, se = se, info = info))
}
