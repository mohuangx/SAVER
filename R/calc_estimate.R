#' Calculate estimate
#'
#' Calculates SAVER estimate
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
#' @param cutoff Maximum absolute correlation to determine whether a gene
#' should be predicted.
#'
#' @param coefs Coefficients of a linear fit of log-squared ratio of
#' largest lambda to lambda of lowest cross-validation error. Used to estimate
#' model with lowest cross-validation error.
#'
#' @param sf Normalized size factor.
#'
#' @param scale.sf Scale of size factor.
#'
#' @param mu Matrix of prior means
#'
#' @param pred.gene.names Names of genes to perform regression prediction.
#'
#' @param pred.cells Index of cells to perform regression prediction.
#'
#' @param null.model Whether to use mean gene expression as prediction.
#'
#' @param nworkers Number of cores registered to parallel backend.
#'
#' @param calc.maxcor Whether to calculate maximum absolute correlation.
#'
#' @param estimates.only Only return SAVER estimates. Default is FALSE.
#'
#' @return A list with the following components
#' \item{\code{est}}{Recovered (normalized) expression}
#' \item{\code{se}}{Standard error of estimates}
#' \item{\code{maxcor}}{Maximum absolute correlation for each gene. 2 if not
#' calculated}
#' \item{\code{lambda.max}}{Smallest value of lambda which gives the null
#' model.}
#' \item{\code{lambda.min}}{Value of lambda from which the prediction model is
#' used}
#' \item{\code{sd.cv}}{Difference in the number of standard deviations in
#' deviance between the model with lowest cross-validation error and the null
#' model}
#' \item{\code{ct}}{Time taken to generate predictions.}
#' \item{\code{vt}}{Time taken to estimate variance.}
#'
#' @rdname calc_estimate
#' @import foreach
#' @export
#'

calc.estimate <- function(x, x.est, cutoff = 0, coefs = NULL, sf, scale.sf,
                          pred.gene.names, pred.cells, null.model, nworkers,
                          calc.maxcor, estimates.only) {
  cs <- min(ceiling(nrow(x)/nworkers), get.chunk(nrow(x), nworkers))
  iterx <- iterators::iter(as.matrix(x), by = "row", chunksize = cs)
  itercount <- iterators::icount(ceiling(iterx$length/iterx$chunksize))
  ix <- NULL; ind <- NULL
  out <- suppressWarnings(
    foreach::foreach(ix = iterx, ind = itercount,
                     .packages = "SAVER", .errorhandling="pass") %dopar% {
      y <- sweep(ix, 2, sf, "/")
      if (calc.maxcor) {
        maxcor <- calc.maxcor(x.est, t(y))
      } else {
        maxcor <- rep(2, nrow(y))
      }

      x.names <- rownames(ix)
      x.est.names <- colnames(x.est)
      est <- matrix(0, nrow(ix), ncol(ix))
      if (!estimates.only) {
        se <- matrix(0, nrow(ix), ncol(ix))
      } else {
        se <- NA
      }
      ct <- rep(0, nrow(ix))
      vt <- rep(0, nrow(ix))
      lambda.max <- rep(0, nrow(ix))
      lambda.min <- rep(0, nrow(ix))
      sd.cv <- rep(0, nrow(ix))

      pred.gene <- (maxcor > cutoff) & (x.names %in% pred.gene.names)
      for (i in 1:nrow(ix)) {
        j <- (ind - 1)*cs + i
        ptc <- Sys.time()
        if (null.model | !pred.gene[i]) {
          pred.out <- list(mean(y[i, pred.cells]), 0, 0, 0)
        } else {
          sameind <- which(x.est.names == x.names[i])
          if (is.null(coefs)) {
            if (length(sameind) == 1) {
              pred.out <- expr.predict(x.est[, -sameind], y[i, ],
                                       pred.cells = pred.cells, seed = j)
            } else {
              pred.out <- expr.predict(x.est, y[i, ],
                                       pred.cells = pred.cells, seed = j)
            }
            lambda.max[i] <- pred.out[[2]]
            lambda.min[i] <- pred.out[[3]]
          } else {
            lambda <- est.lambda(y[i, ], maxcor[i], coefs)
            lambda.max[i] <- lambda[1]
            lambda.min[i] <- lambda[2]
            if (length(sameind) == 1) {
              pred.out <- expr.predict(x.est[, -sameind], y[i, ],
                                       pred.cells = pred.cells,
                                       lambda.max = lambda.max[i],
                                       lambda.min = lambda.min[i])
            } else {
              pred.out <- expr.predict(x.est, y[i, ],
                                       pred.cells = pred.cells,
                                       lambda.max = lambda.max[i],
                                       lambda.min = lambda.min[i])
            }
          }
        }
        ct[i] <- as.numeric(Sys.time()-ptc)
        sd.cv[i] <- pred.out[[4]]
        ptc <- Sys.time()
        post <- calc.post(ix[i, ], pred.out[[1]], sf, scale.sf)
        vt[i] <- as.numeric(Sys.time()-ptc)
        est[i, ] <- post[[1]]
        if (!estimates.only) {
          se[i, ] <- post[[2]]
        }
      }
      list(est, se, maxcor, lambda.max, lambda.min, sd.cv, ct, vt)
    }
  )
  if (length(out[[1]]) != 8) {
    stop(out[[1]])
  }
  est <- do.call(rbind, lapply(out, `[[`, 1))
  if (!estimates.only) {
    se <- do.call(rbind, lapply(out, `[[`, 2))
  } else {
    se <- NA
  }
  maxcor <- unlist(lapply(out, `[[`, 3))
  lambda.max <- unlist(lapply(out, `[[`, 4))
  lambda.min <- unlist(lapply(out, `[[`, 5))
  sd.cv <- unlist(lapply(out, `[[`, 6))
  ct <- unlist(lapply(out, `[[`, 7))
  vt <- unlist(lapply(out, `[[`, 8))
  list(est = est, se = se, maxcor = maxcor, lambda.max = lambda.max,
       lambda.min = lambda.min, sd.cv = sd.cv, ct = ct, vt = vt)
}

#' @rdname calc_estimate
#' @import foreach
#' @export
calc.estimate.mean <- function(x, sf, scale.sf, mu, nworkers, estimates.only) {
  cs <- min(ceiling(nrow(x)/nworkers), get.chunk(nrow(x), nworkers))
  iterx <- iterators::iter(as.matrix(x), by = "row", chunksize = cs)
  itermu <- iterators::iter(mu, by = "row", chunksize = cs)
  itercount <- iterators::icount(ceiling(iterx$length/iterx$chunksize))
  ix <- NULL; ind <- NULL; imu <- NULL
  out <- suppressWarnings(
    foreach::foreach(ix = iterx, imu = itermu, ind = itercount,
                     .packages = "SAVER", .errorhandling="pass") %dopar% {
      y <- sweep(ix, 2, sf, "/")
      maxcor <- rep(0, nrow(y))
      gene.means <- rowMeans(y)
      mu.means <- rowMeans(imu)
      pred <- sweep(imu, 1, rowMeans(y)/rowMeans(imu), "*")
      est <- matrix(0, nrow(ix), ncol(ix))
      if (!estimates.only) {
        se <- matrix(0, nrow(ix), ncol(ix))
      } else {
        se <- NA
      }
      ct <- rep(0, nrow(ix))
      vt <- rep(0, nrow(ix))
      lambda.max <- rep(0, nrow(ix))
      lambda.min <- rep(0, nrow(ix))
      sd.cv <- rep(0, nrow(ix))
      for (i in 1:nrow(ix)) {
        ptc <- Sys.time()
        post <- calc.post(ix[i, ], pred[i, ], sf, scale.sf)
        vt[i] <- as.numeric(Sys.time()-ptc)
        est[i, ] <- post[[1]]
        if (!estimates.only) {
          se[i, ] <- post[[2]]
        }
      }
      list(est, se, maxcor, lambda.max, lambda.min, sd.cv, ct, vt)
    }
  )
  if (length(out[[1]]) != 8) {
    stop(out[[1]])
  }
  est <- do.call(rbind, lapply(out, `[[`, 1))
  if (!estimates.only) {
    se <- do.call(rbind, lapply(out, `[[`, 2))
  } else {
    se <- NA
  }
  maxcor <- unlist(lapply(out, `[[`, 3))
  lambda.max <- unlist(lapply(out, `[[`, 4))
  lambda.min <- unlist(lapply(out, `[[`, 5))
  sd.cv <- unlist(lapply(out, `[[`, 6))
  ct <- unlist(lapply(out, `[[`, 7))
  vt <- unlist(lapply(out, `[[`, 8))
  list(est = est, se = se, maxcor = maxcor, lambda.max = lambda.max,
       lambda.min = lambda.min, sd.cv = sd.cv, ct = ct, vt = vt)
}

#' @rdname calc_estimate
#' @import foreach
#' @export
calc.estimate.null <- function(x, sf, scale.sf, nworkers, estimates.only) {
  cs <- min(ceiling(nrow(x)/nworkers), get.chunk(nrow(x), nworkers))
  iterx <- iterators::iter(as.matrix(x), by = "row", chunksize = cs)
  itercount <- iterators::icount(ceiling(iterx$length/iterx$chunksize))
  ix <- NULL; ind <- NULL
  out <- suppressWarnings(
    foreach::foreach(ix = iterx, ind = itercount, .packages = "SAVER",
                     .errorhandling="pass") %dopar% {
      y <- sweep(ix, 2, sf, "/")
      maxcor <- rep(0, nrow(y))
      pred <- matrix(rowMeans(y), nrow(y), ncol(y))
      est <- matrix(0, nrow(ix), ncol(ix))
      if (!estimates.only) {
        se <- matrix(0, nrow(ix), ncol(ix))
      } else {
        se <- NA
      }
      ct <- rep(0, nrow(ix))
      vt <- rep(0, nrow(ix))
      lambda.max <- rep(0, nrow(ix))
      lambda.min <- rep(0, nrow(ix))
      sd.cv <- rep(0, nrow(ix))
      for (i in 1:nrow(ix)) {
        ptc <- Sys.time()
        post <- calc.post(ix[i, ], pred[i, ], sf, scale.sf)
        vt[i] <- as.numeric(Sys.time()-ptc)
        est[i, ] <- post[[1]]
        if (!estimates.only) {
          se[i, ] <- post[[2]]
        }
      }
      list(est, se, maxcor, lambda.max, lambda.min, sd.cv, ct, vt)
    }
  )
  if (length(out[[1]]) != 8) {
    stop(out[[1]])
  }
  est <- do.call(rbind, lapply(out, `[[`, 1))
  if (!estimates.only) {
    se <- do.call(rbind, lapply(out, `[[`, 2))
  } else {
    se <- NA
  }
  maxcor <- unlist(lapply(out, `[[`, 3))
  lambda.max <- unlist(lapply(out, `[[`, 4))
  lambda.min <- unlist(lapply(out, `[[`, 5))
  sd.cv <- unlist(lapply(out, `[[`, 6))
  ct <- unlist(lapply(out, `[[`, 7))
  vt <- unlist(lapply(out, `[[`, 8))
  list(est = est, se = se, maxcor = maxcor, lambda.max = lambda.max,
       lambda.min = lambda.min, sd.cv = sd.cv, ct = ct, vt = vt)
}
