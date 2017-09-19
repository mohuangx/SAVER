#' Calculates SAVER prediction.
#'
#' Uses \code{cv.glmnet} from the \code{glmnet} package to return the SAVER
#' prediction.
#'
#' The SAVER method starts with predicting the prior mean for each cell for a
#' specific gene. The prediction is performed using the observed normalized
#' gene count as the response and the normalized gene counts of other genes as
#' predictors. \code{cv.glmnet} from the \code{glmnet} package is used to fit
#' the LASSO Poisson regression. The model with the lowest cross-validation
#' error is chosen and the fitted response values are returned and used as the
#' SAVER prediction.
#'
#' @param x A log-normalized expression count matrix of genes to be used in the
#' prediction.
#'
#' @param y A normalized expression count vector of the gene to be predicted.
#'
#' @param pred.cells Index of cells to use for prediction. Default is to use
#' all cells.
#'
#' @param dfmax The number of genes to be included in the prediction. Default
#' is 300.
#'
#' @param nfolds Number of folds to use in the cross-validation. Default is 5.
#'
#' @param nlambda Number of lambda to calculate in cross validation. Default is
#' 5.
#'
#' @param seed Sets the seed for reproducible results.
#'
#'
#' @return A vector of predicted gene expression.
#'
#' @export
expr.predict <- function(x, y, pred.cells = 1:length(y), dfmax = 300,
                         nfolds = 5, nlambda = 50, seed = NULL,
                         verbose = FALSE) {
  if (!is.null(seed))
    set.seed(seed)
  if (sum(y) == 0)
    return(list(mu = rep(0, length(y)), nvar = 0))
  cv <- tryCatch(
    suppressWarnings(glmnet::cv.glmnet(x[pred.cells, ], y[pred.cells],
                                       family="poisson", dfmax = dfmax,
                                       nfolds = nfolds, nlambda = nlambda)),
    error = function(cond) {
      if (verbose)
        message(cond, "\n")
      return(NA)
    }
  )
  if (length(cv) == 1) {
    mu <- rep(mean(y[pred.cells]), length(y))
    nvar <- 0
  } else {
    mu <- c(glmnet::predict.cv.glmnet(cv, newx = x, s = "lambda.min",
                                        type="response"))
    nvar <- cv$nzero[which(cv$lambda == cv$lambda.min)]
  }
  return(list(mu = mu, nvar = nvar))
}
