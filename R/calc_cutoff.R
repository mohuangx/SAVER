
calc.cutoff <- function(x, x.est, sf, npred, pred.cells, nworkers, output.se,
                        verbose, index) {
  cs <- min(ceiling(nrow(x)/nworkers), 10)
  iterx <- iterators::iter(x, by = "row", chunksize = cs)
  itercount <- iterators::icount(ceiling(iterx$length/iterx$chunksize))
  n <- length(pred.cells)
  out <- suppressWarnings(
    foreach::foreach(ix = iterx, ind = itercount,
                     .packages = c("glmnet", "SAVER", "iterators"),
                     .errorhandling="pass") %dopar% {
      y <- sweep(ix, 2, sf, "/")
      if (npred > 100) {
        maxcor <- calc.maxcor(x.est, t(y))
      } else {
        maxcor <- NULL
      }
      
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
        sameind <- which(x.est.names == x.names[i])
        ptc <- proc.time()
        if (length(sameind) == 1) {
          pred.out <- expr.predict(x.est[, -sameind], y[i, ], 
                                   pred.cells = pred.cells,
                                   seed = (ind - 1)*cs + i)
        } else {
          pred.out <- expr.predict(x.est, y[i, ], 
                                   pred.cells = pred.cells,
                                   seed = (ind - 1)*cs + i)
        }
        ct[i] <- (proc.time()-ptc)[3]
        lambda.max[i] <- pred.out[[2]]
        lambda.min[i] <- pred.out[[3]]
        sd.cv[i] <- pred.out[[4]]
        ptc <- proc.time()
        post <- calc.post(ix[i, ], pred.out[[1]], sf, scale.sf)
        vt[i] <- (proc.time()-ptc)[3]
        est[i, ] <- post[[1]]
        if (output.se) {
          se[i, ] <- post[[2]]
        }
      }
      list(est, se, maxcor, lambda.max, lambda.min, sd.cv, ct, vt)
    }
  )
  est <- do.call(rbind, lapply(out, `[[`, 1))
  se <- do.call(rbind, lapply(out, `[[`, 2))
  maxcor <- unlist(lapply(out, `[[`, 3))
  lambda.max <- unlist(lapply(out, `[[`, 4))
  lambda.min <- unlist(lapply(out, `[[`, 5))
  sd.cv <- unlist(lapply(out, `[[`, 6))
  ct <- unlist(lapply(out, `[[`, 7))
  vt <- unlist(lapply(out, `[[`, 8))
  list(est = est, se = se, maxcor = maxcor, lambda.max = lambda.max,
       lambda.min = lambda.min, sd.cv = sd.cv, ct = ct, vt = vt)
}