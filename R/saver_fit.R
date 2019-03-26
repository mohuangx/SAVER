#' Fits SAVER
#'
#' Fits SAVER object
#'
#' The SAVER method starts by estimating the prior mean and variance for the
#' true expression level for each gene and cell. The prior mean is obtained
#' through predictions from a Lasso Poisson regression for each gene
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
#' @param ncores Number of cores to use. Default is 1.
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
#' @param estimates.only Only return SAVER estimates. Default is FALSE.
#'
#' @param mu Matrix of prior means.
#'
#' @return A list with the following components
#' \item{\code{estimate}}{Recovered (normalized) expression}
#' \item{\code{se}}{Standard error of estimates}
#' \item{\code{info}}{Information about fit}
#'
#' @rdname saver_fit
#' @importFrom stats lm
#' @export
#'

saver.fit <- function(x, x.est, do.fast, ncores, sf, scale.sf, pred.genes,
                      pred.cells, null.model, ngenes = nrow(x),
                      ncells = ncol(x), gene.names = rownames(x),
                      cell.names = colnames(x), estimates.only) {
  est <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  if (!estimates.only) {
    se <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  } else {
    se <- NA
  }
  info <- c(list(0), rep(list(rep(0, ngenes)), 6), list(0), list(0), list(0))
  names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min",
                   "sd.cv", "pred.time", "var.time", "cutoff", "lambda.coefs",
                   "total.time")
  info$size.factor <- scale.sf*sf

  nworkers <- ncores
  message("Running SAVER with ", nworkers, " worker(s)")

  pred.genes1 <- pred.genes[Matrix::rowSums(x[pred.genes, , drop = FALSE]) > 0]
  npred1 <- length(pred.genes1)
  npred <- length(pred.genes)
  if (!null.model) {
    message("Calculating predictions for ", npred,
            " genes using ", ncol(x.est), " genes and ", nrow(x.est),
            " cells...")
  } else {
    message("Using means as predictions.")
  }
  st <- Sys.time()
  message("Start time: ", st)
  if (npred1 < ngenes) {
    ind <- c(sample(pred.genes1, npred1), sample((1:ngenes)[-pred.genes1],
                                               ngenes-npred1))
  } else {
    ind <- sample(1:ngenes, ngenes)
  }
  if (do.fast & !null.model) {
    n1 <- min(max(8, nworkers), npred)
    ind1 <- ind[1:n1]
    message("Estimating finish time...")
    t1 <- Sys.time()
    out <- calc.estimate(x[ind1, , drop = FALSE], x.est, cutoff = 0,
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                         pred.cells, null.model, nworkers, calc.maxcor = TRUE,
                         estimates.only)
    est[ind1, ] <- out$est
    if (!estimates.only){
      se[ind1, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind1] <- out[[j+2]]
    }

    n2 <- min(max(100, nworkers), npred)
    t2 <- Sys.time()
    d1 <- mean(info$pred.time[ind[1:n1]] + info$var.time[ind[1:n1]])/nworkers
    perc.pred <- 0.5
    n3 <- min(max(ceiling(n2/perc.pred), nworkers) + n2, npred)
    tdiff <- d1*(n2-n1) + d1*(n3-n2)*perc.pred +
      d1/n1*(npred-n3)*perc.pred/20 +
      (ngenes-npred)*mean(info$var.time[ind[1:n1]])/3/nworkers
    tdiff <- as.difftime(tdiff, units = "secs")
    message("Finished ", n1, "/", ngenes, " genes. Approximate finish time: ",
            Sys.time() + tdiff)


    if (n1 == npred) {
      if (n1 != ngenes) {
        ind6 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind6), " genes.")
        out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                             coefs = NULL, sf, scale.sf,
                             gene.names[pred.genes], pred.cells,
                             null.model = TRUE, nworkers, calc.maxcor = FALSE,
                             estimates.only)
        est[ind6, ] <- out$est
        if (!estimates.only) {
          se[ind6, ] <- out$se
        }
        for (j in 1:6) {
          info[[j+1]][ind6] <- out[[j+2]]
        }
      }
      info[[10]] <- Sys.time() - st
      return(list(estimate = est, se = se, info = info))
    }

    ind2 <- ind[(n1+1):n2]
    message("Calculating max cor cutoff...")
    out <- calc.estimate(x[ind2, , drop = FALSE], x.est, cutoff = 0,
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                         pred.cells, null.model, nworkers, calc.maxcor = TRUE,
                         estimates.only)
    est[ind2, ] <- out$est
    if (!estimates.only) {
      se[ind2, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind2] <- out[[j+2]]
    }

    fit <- lm(sqrt(info$sd.cv[ind[1:n2]]) ~ info$maxcor[ind[1:n2]])
    cutoff <- (0.5 - fit$coefficients[1])/fit$coefficients[2]
    info$cutoff <- unname(cutoff)
    perc.pred <- mean(info$maxcor[ind[1:n2]] > cutoff)

    n3 <- min(max(ceiling(n2/perc.pred), nworkers) + n2, npred)
    ind3 <- ind[(n2+1):n3]

    t3 <- Sys.time()
    d2 <- mean(info$pred.time[ind[1:n2]] + info$var.time[ind[1:n2]])/nworkers
    tdiff <- d2*(n3-n2)*perc.pred + d2*(npred-n3)*perc.pred/20 +
      (ngenes-npred)*mean(info$var.time[ind[1:n2]])/3/nworkers
    tdiff <- as.difftime(tdiff, units = "secs")
    message("Finished ", n2, "/", ngenes, " genes. Approximate finish time: ",
            Sys.time() + tdiff)


    if (n2 == npred) {
      if (n2 != ngenes) {
        ind6 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind6), " genes.")
        out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                             coefs = NULL, sf, scale.sf,
                             gene.names[pred.genes], pred.cells,
                             null.model = TRUE, nworkers, calc.maxcor = FALSE,
                             estimates.only)
        est[ind6, ] <- out$est
        if (!estimates.only) {
          se[ind6, ] <- out$se
        }
        for (j in 1:6) {
          info[[j+1]][ind6] <- out[[j+2]]
        }
      }
      info[[10]] <- Sys.time() - st
      return(list(estimate = est, se = se, info = info))
    }

    message("Calculating lambda coefficients...")
    out <- calc.estimate(x[ind3, , drop = FALSE], x.est, cutoff, coefs = NULL,
                         sf, scale.sf, gene.names[pred.genes], pred.cells,
                         null.model, nworkers, calc.maxcor = TRUE,
                         estimates.only)
    est[ind3, ] <- out$est
    if (!estimates.only) {
      se[ind3, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind3] <- out[[j+2]]
    }

    perc.pred <- mean(info$maxcor[ind[1:n3]] > cutoff)

    t4 <- Sys.time()
    d3 <- mean(info$pred.time[ind[1:n3]] + info$var.time[ind[1:n3]])/nworkers
    tdiff <- d3*(npred-n3)/20 +
      (ngenes-npred)*mean(info$var.time[ind[1:n3]])/3/nworkers
    tdiff <- as.difftime(tdiff, units = "secs")
    message("Finished ", n3, "/", ngenes, " genes. Approximate finish time: ",
            Sys.time() + tdiff)


    if (n3 == npred) {
      if (n3 != ngenes) {
        ind6 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind6), " genes.")
        out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                             coefs = NULL, sf, scale.sf,
                             gene.names[pred.genes], pred.cells,
                             null.model = TRUE, nworkers, calc.maxcor = FALSE,
                             estimates.only)
        est[ind6, ] <- out$est
        if (!estimates.only) {
          se[ind6, ] <- out$se
        }
        for (j in 1:6) {
          info[[j+1]][ind6] <- out[[j+2]]
        }
      }
      info[[10]] <- Sys.time() - st
      return(list(estimate = est, se = se, info = info))
    }

    pred <- which(info$maxcor > cutoff)
    lambda.max <- info$lambda.max[pred]
    lambda.min <- info$lambda.min[pred]
    coefs <- lm(log(lambda.max/lambda.min)^2 ~ info$maxcor[pred])$coefficients
    info[[9]] <- coefs
    cutoff2 <- max(info$cutoff, -coefs[1]/coefs[2])

    n4 <- min(max(ceiling((npred-n3)/4), nworkers) + n3, npred)
    ind4 <- ind[(n3+1):n4]

    message("Predicting remaining genes...")
    out <- calc.estimate(x[ind4, , drop = FALSE], x.est, cutoff2, coefs, sf,
                         scale.sf, gene.names[pred.genes], pred.cells,
                         null.model, nworkers, calc.maxcor = TRUE,
                         estimates.only)

    est[ind4, ] <- out$est
    if (!estimates.only) {
      se[ind4, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind4] <- out[[j+2]]
    }

    t5 <- Sys.time()
    d4 <- difftime(t5, t4, units = "secs")/(n4-n3)
    tdiff <- d4*(npred-n4) +
      (ngenes-npred)*mean(info$var.time[ind[1:n4]])/3/nworkers
    units(tdiff) <- "secs"
    message("Finished ", n4, "/", ngenes, " genes. Approximate finish time: ",
            Sys.time() + tdiff)

    if (n4 == npred) {
      if (n4 != ngenes) {
        ind6 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind6), " genes.")
        out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                             coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                             pred.cells, null.model = TRUE, nworkers,
                             calc.maxcor = FALSE, estimates.only)
        est[ind6, ] <- out$est
        if (!estimates.only) {
          se[ind6, ] <- out$se
        }
        for (j in 1:6) {
          info[[j+1]][ind6] <- out[[j+2]]
        }
      }
      info[[10]] <- Sys.time() - st
      return(list(estimate = est, se = se, info = info))
    }

    n5 <- npred
    ind5 <- ind[(n4+1):n5]

    message("Predicting remaining genes...")
    out <- calc.estimate(x[ind5, , drop = FALSE], x.est, cutoff2, coefs, sf,
                         scale.sf, gene.names[pred.genes], pred.cells,
                         null.model, nworkers, calc.maxcor = TRUE,
                         estimates.only)

    est[ind5, ] <- out$est
    if (!estimates.only) {
      se[ind5, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind5] <- out[[j+2]]
    }

    t6 <- Sys.time()
    tdiff <- (ngenes-npred)*mean(info$var.time[ind[1:n5]])/3/nworkers
    tdiff <- as.difftime(tdiff, units = "secs")

    if (n5 != ngenes) {
      message("Finished ", n5, "/", ngenes, " genes. Approximate finish time: ",
              Sys.time() + tdiff)
      ind6 <- ind[(npred+1):ngenes]
      message("Estimating remaining ", length(ind6), " genes.")
      out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                           coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                           pred.cells, null.model = TRUE, nworkers,
                           calc.maxcor = FALSE, estimates.only)
      est[ind6, ] <- out$est
      if (!estimates.only) {
        se[ind6, ] <- out$se
      }
      for (j in 1:6) {
        info[[j+1]][ind6] <- out[[j+2]]
      }
    }
    info[[10]] <- Sys.time() - st
    return(list(estimate = est, se = se, info = info))
  } else {
    n1 <- min(max(8, nworkers), npred)
    ind1 <- ind[1:n1]
    message("Estimating finish time...")
    t1 <- Sys.time()
    out <- calc.estimate(x[ind1, , drop = FALSE], x.est, cutoff = 0,
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                         pred.cells, null.model, nworkers, calc.maxcor = FALSE,
                         estimates.only)
    est[ind1, ] <- out$est
    if (!estimates.only) {
      se[ind1, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind1] <- out[[j+2]]
    }

    t2 <- Sys.time()
    d1 <- mean(info$pred.time[ind[1:n1]] + info$var.time[ind[1:n1]])/nworkers
    tdiff <- d1*(npred-n1) +
      (ngenes-npred)*mean(info$var.time[ind[1:n1]])/3/nworkers
    tdiff <- as.difftime(tdiff, units = "secs")
    message("Finished ", n1, "/", ngenes, " genes. Approximate finish time: ",
            Sys.time() + tdiff)

    if (n1 == npred) {
      if (n1 != ngenes) {
        ind6 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind6), " genes.")
        out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                             coefs = NULL, sf, scale.sf,
                             gene.names[pred.genes], pred.cells,
                             null.model = TRUE, nworkers, calc.maxcor = FALSE,
                             estimates.only)
        est[ind6, ] <- out$est
        if (!estimates.only) {
          se[ind6, ] <- out$se
        }
        for (j in 1:6) {
          info[[j+1]][ind6] <- out[[j+2]]
        }
      }
      info[[10]] <- Sys.time() - st
      return(list(estimate = est, se = se, info = info))
    }

    n2 <- min(ceiling((npred-n1)/4) + n1, npred)
    ind2 <- ind[(n1+1):n2]
    message("Predicting remaining genes...")
    out <- calc.estimate(x[ind2, , drop = FALSE], x.est, cutoff = 0,
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                         pred.cells, null.model, nworkers, calc.maxcor = FALSE,
                         estimates.only)
    est[ind2, ] <- out$est
    if (!estimates.only) {
      se[ind2, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind2] <- out[[j+2]]
    }

    t3 <- Sys.time()
    d2 <- difftime(t3, t1, units = "secs")/n2
    tdiff <- d2*(npred-n2) +
      (ngenes-npred)*mean(info$var.time[ind[1:n2]])/3/nworkers
    message("Finished ", n2, "/", ngenes, " genes. Approximate finish time: ",
            Sys.time() + tdiff)

    if (n2 == npred) {
      if (n2 != ngenes) {
        ind6 <- ind[(npred+1):ngenes]
        message("Estimating remaining ", length(ind6), " genes.")
        out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                             coefs = NULL, sf, scale.sf,
                             gene.names[pred.genes], pred.cells,
                             null.model = TRUE, nworkers, calc.maxcor = FALSE,
                             estimates.only)
        est[ind6, ] <- out$est
        if (!estimates.only) {
          se[ind6, ] <- out$se
        }
        for (j in 1:6) {
          info[[j+1]][ind6] <- out[[j+2]]
        }
      }
      info[[10]] <- Sys.time() - st
      return(list(estimate = est, se = se, info = info))
    }

    n3 <- npred
    ind3 <- ind[(n2+1):n3]

    message("Predicting remaining genes...")
    out <- calc.estimate(x[ind3, , drop = FALSE], x.est, cutoff = 0,
                         coefs = NULL, sf, scale.sf, gene.names[pred.genes],
                         pred.cells, null.model, nworkers, calc.maxcor = FALSE,
                         estimates.only)
    est[ind3, ] <- out$est
    if (!estimates.only) {
      se[ind3, ] <- out$se
    }
    for (j in 1:6) {
      info[[j+1]][ind3] <- out[[j+2]]
    }

    t4 <- Sys.time()
    tdiff <- (ngenes-npred)*mean(info$var.time[ind[1:n3]])/3/nworkers
    tdiff <- as.difftime(tdiff, units = "secs")

    if (n3 != ngenes) {
      message("Finished ", n3, "/", ngenes, " genes. Approximate finish time: ",
              Sys.time() + tdiff)
      ind6 <- ind[(npred+1):ngenes]
      message("Estimating remaining ", length(ind6), " genes.")
      out <- calc.estimate(x[ind6, , drop = FALSE], x.est, cutoff = 0,
                           coefs = NULL, sf, scale.sf,
                           gene.names[pred.genes], pred.cells,
                           null.model = TRUE, nworkers, calc.maxcor = FALSE,
                           estimates.only)
      est[ind6, ] <- out$est
      if (!estimates.only) {
        se[ind6, ] <- out$se
      }
      for (j in 1:6) {
        info[[j+1]][ind6] <- out[[j+2]]
      }
    }
    info[[10]] <- Sys.time() - st
    return(list(estimate = est, se = se, info = info))
  }
}

#' @rdname saver_fit
#' @export
saver.fit.mean <- function(x, ncores, sf, scale.sf, mu, ngenes = nrow(x),
                           ncells = ncol(x), gene.names = rownames(x),
                           cell.names = colnames(x), estimates.only) {
  out <- list()
  out$estimate <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  if (!estimates.only) {
    out$se <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  } else {
    out$se <- NA
  }
  out$info <- c(list(0), rep(list(rep(0, ngenes)), 6), list(0), list(0), list(0))
  names(out$info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min",
                   "sd.cv", "pred.time", "var.time", "cutoff", "lambda.coefs",
                   "total.time")
  out$info$size.factor <- scale.sf*sf

  nworkers <- ncores
  message("Running SAVER given prior means with ", nworkers, " worker(s)")

  st <- Sys.time()
  message("Start time: ", st)
  ind <- sample(1:ngenes, ngenes)

  n1 <- min(max(100, nworkers), ngenes)
  ind1 <- ind[1:n1]
  message("Estimating finish time...")
  t1 <- Sys.time()
  results <- calc.estimate.mean(x[ind1, , drop = FALSE], sf, scale.sf,
                            mu[ind1, , drop = FALSE], nworkers, estimates.only)
  out$estimate[ind1, ] <- results$est
  if (!estimates.only) {
    out$se[ind1, ] <- results$se
  }
  for (j in 1:6) {
    out$info[[j+1]][ind1] <- results[[j+2]]
  }

  if (n1 == ngenes) {
    out$info[[10]] <- Sys.time() - st
    return(out)
  }
  
  d1 <- mean(out$info$var.time[ind[1:n1]])/nworkers
  tdiff <- d1*(ngenes-n1)
  tdiff <- as.difftime(tdiff, units = "secs")
  message("Finished ", n1, "/", ngenes, " genes. Approximate finish time: ",
          Sys.time() + tdiff)
  
  ngenes.left <- ngenes-n1
  total.elem <- ngenes.left*ncells
  
  nsplit <- max(4, ceiling(total.elem/(2^31-1)))
  
  split.ind <- ceiling(seq(n1, ngenes, length.out = nsplit+1))
  
  t1 <- Sys.time()
  for (i in 1:(length(split.ind)-1)) {
    out <- update.output(calc.estimate.mean, ind, split.ind[i]+1, 
                         split.ind[i+1], out, x, sf, scale.sf, mu, 
                         nworkers, estimates.only)
    t2 <- Sys.time()
    d1 <- difftime(t2, t1, units = "secs")/(split.ind[i+1]-n1)
    tdiff <- d1*(ngenes-split.ind[i+1])
    tdiff <- as.difftime(tdiff, units = "secs")
    message("Finished ", split.ind[i+1], "/", ngenes, 
            " genes. Approximate finish time: ",
            Sys.time() + tdiff)
  }

  out$info[[10]] <- Sys.time() - st
  return(out)
}

#' @rdname saver_fit
#' @export
saver.fit.null <- function(x, ncores, sf, scale.sf, ngenes = nrow(x),
                           ncells = ncol(x), gene.names = rownames(x),
                           cell.names = colnames(x), estimates.only) {
  est <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  if (!estimates.only) {
    se <- matrix(0, ngenes, ncells, dimnames = list(gene.names, cell.names))
  } else {
    se <- NA
  }
  info <- c(list(0), rep(list(rep(0, ngenes)), 6), list(0), list(0), list(0))
  names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min",
                   "sd.cv", "pred.time", "var.time", "cutoff", "lambda.coefs",
                   "total.time")
  info$size.factor <- scale.sf*sf

  nworkers <- ncores
  message("Running SAVER given null model with ", nworkers, " worker(s)")

  st <- Sys.time()
  message("Start time: ", st)
  ind <- sample(1:ngenes, ngenes)

  n1 <- min(max(100, nworkers), ngenes)
  ind1 <- ind[1:n1]
  message("Estimating finish time...")
  t1 <- Sys.time()
  out <- calc.estimate.null(x[ind1, , drop = FALSE], sf, scale.sf, nworkers,
                            estimates.only)
  est[ind1, ] <- out$est
  if (!estimates.only) {
    se[ind1, ] <- out$se
  }
  for (j in 1:6) {
    info[[j+1]][ind1] <- out[[j+2]]
  }

  if (n1 == ngenes) {
    info[[10]] <- Sys.time() - st
    return(list(estimate = est, se = se, info = info))
  }

  n2 <- min(ceiling((ngenes-n1)/4 + n1), ngenes)
  t2 <- Sys.time()
  d1 <- mean(info$var.time[ind[1:n1]])/nworkers
  tdiff <- d1*(ngenes-n1)
  tdiff <- as.difftime(tdiff, units = "secs")
  message("Finished ", n1, "/", ngenes, " genes. Approximate finish time: ",
          Sys.time() + tdiff)

  ind2 <- ind[(n1+1):n2]
  message("Estimating remaining genes...")
  out <- calc.estimate.null(x[ind2, , drop = FALSE], sf, scale.sf, nworkers,
                            estimates.only)
  est[ind2, ] <- out$est
  if (!estimates.only) {
    se[ind2, ] <- out$se
  }
  for (j in 1:6) {
    info[[j+1]][ind2] <- out[[j+2]]
  }

  if (n2 == ngenes) {
    info[[10]] <- Sys.time() - st
    return(list(estimate = est, se = se, info = info))
  }

  n3 <- ngenes
  t3 <- Sys.time()
  d2 <- difftime(t3, t1, units = "secs")/n2
  tdiff <- d2*(n3-n2)
  message("Finished ", n2, "/", ngenes, " genes. Approximate finish time: ",
          Sys.time() + tdiff)

  ind3 <- ind[(n2+1):n3]
  message("Estimating remaining genes...")
  out <- calc.estimate.null(x[ind3, , drop = FALSE], sf, scale.sf, nworkers,
                            estimates.only)
  est[ind3, ] <- out$est
  if (!estimates.only) {
    se[ind3, ] <- out$se
  }
  for (j in 1:6) {
    info[[j+1]][ind3] <- out[[j+2]]
  }
  info[[10]] <- Sys.time() - st
  return(list(estimate = est, se = se, info = info))
}


