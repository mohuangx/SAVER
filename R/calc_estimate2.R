
calc.estimate <- function(x, x.est, cutoff = 0, coefs = NULL, sf, scale.sf, 
                          pred.gene.names, pred.cells, null.model, nworkers, 
                          verbose) {
  cs <- min(ceiling(nrow(x)/nworkers), 100)
  iterx <- iterators::iter(x, by = "row", chunksize = cs)
  itercount <- iterators::icount(ceiling(iterx$length/iterx$chunksize))
  out <- suppressWarnings(
    foreach::foreach(ix = iterx, ind = itercount,
                     .packages = c("glmnet", "SAVER", "iterators"),
                     .errorhandling="pass") %dopar% {
      y <- sweep(ix, 2, sf, "/")
      if (length(pred.gene.names) > 100) {
        maxcor <- SAVER::calc.maxcor(x.est, t(y))
      } else {
        maxcor <- rep(2, nrow(y))
      }
      
      x.names <- rownames(ix)
      x.est.names <- colnames(x.est)
      est <- matrix(0, nrow(ix), ncol(ix))
      se <- matrix(0, nrow(ix), ncol(ix))
      ct <- rep(0, nrow(ix))
      vt <- rep(0, nrow(ix))
      lambda.max <- rep(0, nrow(ix))
      lambda.min <- rep(0, nrow(ix))
      sd.cv <- rep(0, nrow(ix))
      
      pred.gene <- (maxcor > cutoff) & (x.names %in% pred.gene.names)
      for (i in 1:nrow(ix)) {
        ptc <- proc.time()
        if (null.model | !pred.gene[i]) {
          pred.out <- list(mean(y[i, pred.cells]), 0, 0, 0)
        } else {
          sameind <- which(x.est.names == x.names[i])
          if (is.null(coefs)) {
            if (length(sameind) == 1) {
              pred.out <- SAVER::expr.predict(x.est[, -sameind], y[i, ], 
                                              pred.cells = pred.cells,
                                              seed = (ind - 1)*cs + i)
            } else {
              pred.out <- SAVER::expr.predict(x.est, y[i, ], 
                                              pred.cells = pred.cells,
                                              seed = (ind - 1)*cs + i)
            }
          } else {
            lambda <- est.lambda(y[i, ], maxcor[i], coefs)
            lambda.max[i] <- lambda[1]
            lambda.min[i] <- lambda[2]
            if (length(sameind) == 1) {
              pred.out <- SAVER::expr.predict(x.est[, -sameind], y[i, ], 
                                              pred.cells = pred.cells,
                                              lambda.max = lambda.max[i],
                                              lambda.min = lambda.min[i])
            } else {
              pred.out <- SAVER::expr.predict(x.est, y[i, ], 
                                              pred.cells = pred.cells,
                                              lambda.max = lambda.max[i],
                                              lambda.min = lambda.min[i])
            }
          }
        }
        ct[i] <- (proc.time()-ptc)[3]
        sd.cv[i] <- pred.out[[4]]
        ptc <- proc.time()
        post <- SAVER::calc.post(ix[i, ], pred.out[[1]], sf, scale.sf)
        vt[i] <- (proc.time()-ptc)[3]
        est[i, ] <- post[[1]]
        se[i, ] <- post[[2]]
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