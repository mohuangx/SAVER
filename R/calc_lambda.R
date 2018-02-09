
calc.lambda <- function(x, x.est, npred, pred.cells, nworkers, output.se,
                        verbose, index) {
  cs <- min(ceiling(nrow(x)/nworkers), 10)
  iterx <- iterators::iter(x, by = "row", chunksize = cs)
  itercount <- iterators::icount(ceiling(iterx$length/iterx$chunksize))
  n <- length(pred.cells)
  out <- suppressWarnings(
    foreach::foreach(ix = iterx, ind = itercount,
                     .packages = c("glmnet", "SAVER", "iterators"),
                     .errorhandling="pass") %dopar% {
      maxcor <- calc.maxcor(x.est, t(ix))
      x.names <- rownames(ix)
      x.est.names <- colnames(x.est)
      est <- matrix(0, nrow(ix), ncol(ix))
      if (output.se) {
        se <- matrix(0, nrow(ix), ncol(ix))
      } else {
        se <- NULL
      }
      lambda.min <- rep(0, nrow(ix))
      sd.cv <- rep(0, nrow(ix))
      
      for (i in 1:nrow(ix)) {
        y <- ix[i, pred.cells]/sf[pred.cells]
        if (maxcor[i] > cutoff) {
          sameind <- which(x.est.names == x.names[i])
          if (length(sameind) == 1) {
            pred.out <- expr.predict.cv(x.est[pred.cells, -sameind], y,
                                        seed = (ind - 1)*cs + i)
          } else {
            pred.out <- expr.predict.cv(x.est[pred.cells, ], y,
                                        seed = (ind - 1)*cs + i)
          }
        } else {
          pred.out <- list(mean(y), 0, 0)
        }

        lambda.min[i] <- pred.out[[2]]
        sd.cv[i] <- pred.out[[3]]
        post <- calc.post(ix[i, ], pred.out[[1]], sf, scale.sf)
        est[i, ] <- post[[1]]
        if (output.se) {
          se[i, ] <- post[[2]]
        }
      }
      list(est, se, maxcor, lambda.min, sd.cv)
    }
  )
  est <- do.call(rbind, lapply(out, `[[`, 1))
  se <- do.call(rbind, lapply(out, `[[`, 2))
  maxcor <- unlist(lapply(out, `[[`, 3))
  lambda.min <- unlist(lapply(out, `[[`, 4))
  sd.cv <- unlist(lapply(out, `[[`, 5))
  fit <- lm(sqrt(sd.cv) ~ maxcor)
  cutoff <- (0.5 - fit$coefficients[1])/fit$coefficients[2]
  list(est = est, se = se, maxcor = maxcor, lambda.min = lambda.min, 
       sd.cv = sd.cv, cutoff = cutoff)
}