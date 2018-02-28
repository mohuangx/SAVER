
clean.data <- function(x) {
  if (!is.matrix(x)) {
    x <- as.matrix(x)
    message("Converting x to matrix.")
    if (!is.numeric(x)) {
      stop("Make sure x is numeric.")
    }
  }
  np <- dim(x)
  if (is.null(np) | (np[2] <= 1))
    stop("x should be a matrix with 2 or more columns")
  if (min(colSums(x)) == 0) {
    nzerocells <- sum(colSums(x) == 0)
    x <- x[, colSums(x) != 0]
    message("Removing ", nzerocells, " cell(s) with zero expression.")
  }
  if (is.null(rownames(x))) {
    rownames(x) <- 1:np[1]
  }
  x
}

calc.size.factor <- function(x, size.factor, ncells) {
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
    pred.genes <- 1:ngenes
  } else if (npred < ngenes) {
    pred.genes <- order(rowMeans(x), decreasing = TRUE)[1:npred]
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