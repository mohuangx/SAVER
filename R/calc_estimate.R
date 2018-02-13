
calc.estimate <- function(x, x.est, cutoff, fit, sf, npred, pred.cells, 
                          nworkers, output.se, verbose, index) {
  cs <- min(ceiling(nrow(x)/nworkers), 10)
  iterx <- iterators::iter(x, by = "row", chunksize = cs)
  itercount <- iterators::icount(ceiling(iterx$length/iterx$chunksize))
  n <- length(pred.cells)
  a1 <- system.time(out <- suppressWarnings(
    foreach::foreach(ix = iterx, ind = itercount,
                     .packages = c("glmnet", "SAVER", "iterators"),
                     .errorhandling="pass") %dopar% {
      ptc <- proc.time()
      maxcor <- calc.maxcor(x.est, t(sweep(ix, 2, sf, "/")))
      x.names <- rownames(ix)
      x.est.names <- colnames(x.est)
      est <- matrix(0, nrow(ix), ncol(ix))
      if (output.se) {
        se <- matrix(0, nrow(ix), ncol(ix))
      } else {
        se <- NULL
      }
      ct <- rep(0, nrow(ix))
      vt <- rep(0, nrow(ix))
      lambda.max <- rep(0, nrow(ix))
      lambda.min <- rep(0, nrow(ix))
      sd.cv <- rep(0, nrow(ix))
      
      for (i in 1:nrow(ix)) {
        y <- ix[i, pred.cells]/sf[pred.cells]
        if (maxcor[i] > cutoff) {
          sameind <- which(x.est.names == x.names[i])
          lambda <- est.lambda(y, maxcor[i], fit$coefficients)
          lambda.max[i] <- lambda[1]
          lambda.min[i] <- lambda[2]
          if (length(sameind) == 1) {
            ct[i] <- system.time(pred.out <- expr.predict(
              x.est[pred.cells, -sameind], y, lambda.max = lambda.max[i],
              lambda.min = lambda.min[i]))[3]
          } else {
            ct[i] <- system.time(pred.out <- expr.predict(
              x.est[pred.cells, ], y, lambda.max = lambda.max[i],
              lambda.min = lambda.min[i]))[3]
          }
        } else {
          pred.out <- list(mean(y[pred.cells]), 0, 0, 0)
        }
        sd.cv[i] <- pred.out[[4]]
        vt[i] <- system.time(post <- calc.post(ix[i, ], pred.out[[1]], sf, 
                                               scale.sf))[3]
        est[i, ] <- post[[1]]
        if (output.se) {
          se[i, ] <- post[[2]]
        }
      }
      a3 <- (proc.time()-ptc)[3]
      list(est, se, maxcor, lambda.max, lambda.min, sd.cv, ct, vt, a3)
    }
  ))[3]
  ptc <- proc.time()
  est <- do.call(rbind, lapply(out, `[[`, 1))
  se <- do.call(rbind, lapply(out, `[[`, 2))
  maxcor <- unlist(lapply(out, `[[`, 3))
  lambda.max <- unlist(lapply(out, `[[`, 4))
  lambda.min <- unlist(lapply(out, `[[`, 5))
  sd.cv <- unlist(lapply(out, `[[`, 6))
  ct <- unlist(lapply(out, `[[`, 7))
  vt <- unlist(lapply(out, `[[`, 8))
  a2 <- proc.time()-ptc
  list(est = est, se = se, maxcor = maxcor, lambda.max = lambda.max,
       lambda.min = lambda.min, sd.cv = sd.cv, ct = ct, vt = vt,
       a1 = a1, a2 = a2, a3 = a3)
}