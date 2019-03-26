
clean.data <- function(x) {
  if (!(grepl("matrix", class(x), ignore.case = TRUE))) {
    x <- Matrix::Matrix(as.matrix(x))
    message("Converting x to matrix.")
    if (!is.numeric(x)) {
      warning("Make sure x is numeric.")
    }
  }
  np <- dim(x)
  size <- as.numeric(np[1])*as.numeric(np[2])
  if(size > 2^31-1){
    inds <- split(1:np[2], ceiling(1:np[2]/1000))
    for(i in 1:length(inds)){
      x[, inds[[i]]][x[, inds[[i]]] < 0.001] <- 0
    }
  } else {
    x[x < 0.001] <- 0
  }
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  if (min(Matrix::colSums(x)) == 0) {
    nzerocells <- sum(Matrix::colSums(x) == 0)
    x <- x[, Matrix::colSums(x) != 0]
    message("Removing ", nzerocells, " cell(s) with zero expression.")
  }
  if (is.null(rownames(x))) {
    rownames(x) <- 1:np[1]
  }
  x
}

check.mu <- function(x, mu) {
  if (!is.matrix(mu)) {
    mu <- as.matrix(mu)
    message("Converting mu to matrix.")
    if (!is.numeric(mu)) {
      stop("Make sure mu is numeric.")
    }
  }
  np <- dim(x)
  npmu <- dim(mu)
  if (sum(np == npmu) != 2) {
    stop("x and mu must have same dimensions")
  }
  if (min(Matrix::colSums(x)) == 0) {
    nzerocells <- sum(Matrix::colSums(x) == 0)
    mu <- mu[, Matrix::colSums(x) != 0]
  }
  mu
}

update.output <- function(f, ind, start, stop, out, x, sf, scale.sf, mu, 
                          nworkers, estimates.only) {
  n <- stop-start+1
  ind1 <- ind[start:stop]
  results <- f(x[ind1, , drop = FALSE], sf, scale.sf,
               mu[ind1, , drop = FALSE], nworkers, estimates.only)
  out$estimate[ind1, ] <- results$est
  if (!estimates.only) {
    out$se[ind1, ] <- results$se
  }
  for (j in 1:6) {
    out$info[[j+1]][ind1] <- results[[j+2]]
  }
  return(out)
}


calc.size.factor <- function(x, size.factor, ncells) {
  if (is.null(size.factor)) {
    sf <- Matrix::colSums(x)/mean(Matrix::colSums(x))
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
  list(unname(sf), scale.sf)
}

get.pred.cells <- function(pred.cells, ncells) {
  if (!is.null(pred.cells)) {
    if (min(pred.cells) < 1 |
        max(pred.cells) > ncells) {
      stop("pred.cells must be column indices of x")
    }
  } else {
    pred.cells <- 1:ncells
  }
  pred.cells
}

get.pred.genes <- function(x, pred.genes, npred, ngenes) {
  if (!is.null(pred.genes)) {
    if (min(pred.genes) < 1 |
        max(pred.genes) > ngenes) {
      stop("pred.genes must be row indices of x")
    }
  } else if (is.null(npred)) {
    pred.genes <- (1:ngenes)[Matrix::rowSums(x) != 0]
  } else if (npred < ngenes) {
    npred <- min(sum(Matrix::rowSums(x) != 0), npred)
    pred.genes <- order(Matrix::rowMeans(x), decreasing = TRUE)[1:npred]
  } else {
    stop("npred must be less than number of rows in x")
  }
  pred.genes
}

# split genes for parallel
get.chunk <- function(len, n) {
  cs <- c(rbind(100:150, 99:49))
  r <- rep(0, length(cs))
  for (i in 1:length(cs)) {
    r[i] <- ceiling(len/cs[i]) %% n
    if (r[i] == 0) return(cs[i])
  }
  return(cs[which.max(r)])
}

# estimate lambda.min
est.lambda <- function(x, r, coefs) {
  n <- length(x)
  lambda.max <- r*sd(x)*sqrt((n-1)/n)
  lambda.min <- lambda.max*exp(-sqrt(max(0, coefs[1]+coefs[2]*r)))
  c(lambda.max, lambda.min)
}

# old version of sample.saver
sample.saver.old <- function(x, rep = 1, efficiency.known = FALSE,
                             seed = NULL) {
  ncells <- ncol(x$estimate)
  ngenes <- nrow(x$estimate)
  cell.names <- colnames(x$estimate)
  gene.names <- rownames(x$estimate)
  rep <- as.integer(rep)
  if (rep <= 0) {
    stop("rep must be a positive integer.")
  }
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (rep == 1) {
    if (efficiency.known) {
      samp <- t(sapply(1:ngenes, function(i)
        rnbinom(ncells, mu = x$estimate[i, ], size = x$alpha[i, ])))
    } else {
      samp <- t(sapply(1:ngenes, function(i)
        rgamma(ncells, x$alpha[i, ], x$beta[i, ])))
      samp <- round(samp, 3)
    }
    rownames(samp) <- gene.names
    colnames(samp) <- cell.names
  } else {
    samp <- vector("list", rep)
    for (j in 1:rep) {
      if (efficiency.known) {
        samp[[j]] <- t(sapply(1:ngenes, function(i)
          rnbinom(ncells, mu = x$estimate[i, ], size = x$alpha[i, ])))
      } else {
        samp[[j]] <- t(sapply(1:ngenes, function(i)
          rgamma(ncells, x$alpha[i, ], x$beta[i, ])))
        samp[[j]] <- round(samp[[j]], 3)
      }
      rownames(samp[[j]]) <- gene.names
      colnames(samp[[j]]) <- cell.names
    }
  }
  return(samp)
}

# old version of cor.genes and cor.cells
cor.genes.old <- function(x, cor.mat = NULL) {
  if (is.null(cor.mat)) {
    message("Calculating correlation matrix...")
    cor.mat <- cor(t(x$estimate))
  }
  ngenes <- nrow(x$estimate)
  adj.vec <- rep(0, ngenes)
  for (i in 1:ngenes) {
    adj.vec[i] <- sqrt(var(x$estimate[i, ], na.rm = TRUE)/
                         (var(x$estimate[i, ], na.rm = TRUE) +
                            mean(x$alpha[i, ]/x$beta[i, ]^2, na.rm = TRUE)))
  }
  adj.mat <- outer(adj.vec, adj.vec)
  cor.adj <- adj.mat*cor.mat
  return(cor.adj)
}

cor.cells.old <- function(x, cor.mat = NULL) {
  if (is.null(cor.mat)) {
    message("Calculating correlation matrix...")
    cor.mat <- cor(x$estimate)
  }
  ncells <- ncol(x$estimate)
  adj.vec <- rep(0, ncells)
  for (i in 1:ncells) {
    adj.vec[i] <- sqrt(var(x$estimate[, i], na.rm = TRUE)/
                         (var(x$estimate[, i], na.rm = TRUE) +
                            mean(x$alpha[, i]/x$beta[, i]^2, na.rm = TRUE)))
  }
  adj.mat <- outer(adj.vec, adj.vec)
  cor.adj <- adj.mat*cor.mat
  return(cor.adj)
}

# old version of combine.saver
combine.saver.old <- function(saver.list) {
  est <- do.call(rbind, lapply(saver.list, `[[`, 1))
  alpha <- do.call(rbind, lapply(saver.list, `[[`, 2))
  beta <- do.call(rbind, lapply(saver.list, `[[`, 3))
  info <- vector("list", length(saver.list[[1]]$info))
  names(info) <- names(saver.list[[1]]$info)
  info[[1]] <- saver.list[[1]]$info$size.factor
  info.list <- lapply(saver.list, `[[`, 4)
  for (i in 2:3) {
    info[[i]] <- do.call(c, lapply(info.list, `[[`, i))
  }
  for (i in 4:length(info)) {
    info[[i]] <- saver.list[[1]]$info[[i]]
  }
  out <- list(estimate = est, alpha = alpha, beta = beta, info = info)
  class(out) <- "saver"
  out
}

