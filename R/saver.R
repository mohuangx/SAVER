

saver <- function(x, size.factor = NULL, npred = NULL, pred.genes = NULL,
                  pred.genes.only = FALSE, parallel = FALSE) {
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  if (min(rowSums(x)) == 0 | min(colSums(x)) == 0)
    stop("Please remove genes or cells with all zero expression")
  ngenes <- as.integer(np[1])
  ncells <- as.integer(np[2])
  if (is.null(size.factor)) {
    sf <- colSums(x)/mean(colSums(x))
    scale.sf <- 1
  } else if (length(size.factor) == ncells) {
    sf <- size.factor/mean(size.factor)
    scale.sf <- mean(size.factor)
  } else if (size.factor == 1) {
    sf <- rep(1, ncells)
    scale.sf <- 1
  } else if (min(size.factor) <= 0) {
    stop("Size factor must be greater than 0")
  } else {
    stop("Not a valid size factor")
  }
  if (!is.null(pred.genes)) {
    if (!is.integer(pred.genes) | min(pred.genes) < 1 |
        max(pred.genes) > as.integer(np[1])) {
      stop("pred.genes must be row indices of x")
    }
    npred <- length(pred.genes)
  } else if (is.null(npred)) {
    npred <- ngenes
    pred.genes <- 1:ngenes
  } else if (npred < ngenes) {
    pred.genes <- order(rowMeans(x), decreasing = TRUE)[1:npred]
  } else {
    stop("npred must be less than number of rows in x")
  }
  if (pred.genes.only)
    ngenes <- npred
  x.norm <- sweep(x, 2, sf, "/")
  x.est <- log(sweep(x.norm, 2, 1/sf, "+"))
  gene.means <- rowMeans(x.norm)
  message("calculating predictions...")
  t1 <- Sys.time()
  cvt <- expr.predict(t(x.est[-pred.genes[1], ]), x[pred.genes[1], ]/sf)
  t2 <- Sys.time()
  if (parallel) {
    nworkers <- foreach::getDoParWorkers()
    t3 <- (t2-t1)*npred/nworkers + 0.05*ngenes
    units(t3) <- "mins"
    message("Expected time to finish: ", round(t3))
    if (nworkers == 1) {
      message("Only one worker assigned! Running sequentially...")
      mu.par <- matrix(0, npred, ncells)
      for (i in 1:length(pred.genes)) {
        mu.par[i, ] <- expr.predict(t(x.est[-pred.genes[i], ]),
                                    x[pred.genes[i], ]/sf,
                                    seed = pred.genes[i])
        message(i)
      }
    } else {
      message("Running in parallel: ", nworkers, " workers")
      gene.list <- chunk2(pred.genes, nworkers)
      mu.par <- foreach::foreach(i = 1:nworkers, .combine = rbind,
                                 .packages = c("glmnet", "SAVER")) %dopar% {
        mu.temp <- matrix(0, length(gene.list[[i]]), ncells)
        for (j in 1:length(gene.list[[i]])) {
          mu.temp[j, ] <- expr.predict(t(x.est[-gene.list[[i]][j], ]),
                                       x[gene.list[[i]][j], ]/sf,
                                       seed = gene.list[[i]][j])
        }
        return(mu.temp)
      }
    }

  } else {
    t3 <- (t2-t1)*npred + 0.05*ngenes
    units(t3) <- "mins"
    message("Expected time to finish: ", round(t3))
    mu.par <- matrix(0, npred, ncells)
    for (i in 1:length(pred.genes)) {
      mu.par[i, ] <- expr.predict(t(x.est[-pred.genes[i], ]),
                                  x[pred.genes[i], ]/sf, seed = pred.genes[i])
      message(i)
    }
  }
  if (pred.genes.only) {
    out <- lapply(1:4, function(x) matrix(0, npred, ncells))
    mu <- mu.par
    for (i in 1:length(pred.genes)) {
      post <- calc.post(x[pred.genes[i], ], mu[i, ], sf)
      out[[1]][i, ] <- post$estimate
      out[[2]][i, ] <- post$alpha
      out[[3]][i, ] <- post$beta
    }
    out[[4]] <- mu
    gene.names <- rownames(x)[pred.genes]
    cell.names <- colnames(x)
  } else {
    out <- lapply(1:4, function(x) matrix(0, ngenes, ncells))
    mu <- matrix(gene.means, ngenes, ncells)
    mu[pred.genes, ] <- mu.par
    for (i in 1:ngenes) {
      post <- calc.post(x[i, ], mu[i, ], sf)
      out[[1]][i, ] <- post$estimate
      out[[2]][i, ] <- post$alpha
      out[[3]][i, ] <- post$beta
    }
    out[[4]] <- mu
    gene.names <- rownames(x)
    cell.names <- colnames(x)
  }
  out.named <- lapply(out, function(x) {rownames(x) <- gene.names;
                                        colnames(x) <- cell.names; x})
  names(out.named) <- c("estimate", "alpha", "beta", "predicted")
  return(out.named)
}

chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
