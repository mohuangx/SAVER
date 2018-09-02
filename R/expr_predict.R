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
#' @param seed Sets the seed for reproducible results.
#'
#' @param lambda.max Maximum value of lambda which gives null model.
#'
#' @param lambda.min Value of lambda from which the prediction model is
#' used
#'
#' @return A vector of predicted gene expression.
#'
#' @importFrom stats sd
#' @export
expr.predict <- function(x, y, pred.cells = 1:length(y), seed = NULL,
                         lambda.max = NULL, lambda.min = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  if (sd(y) == 0)
    return(list(rep(mean(y[pred.cells]), length(y)), 0, 0, 0))
  if (is.null(lambda.max)) {
    cv <- tryCatch(
      suppressWarnings(glmnet::cv.glmnet(x[pred.cells, ], y[pred.cells],
                                         family="poisson", dfmax = 300,
                                         nfolds = 5)),
      error = function(cond) {
        return(NA)
      }
    )
    if (length(cv) == 1) {
      mu <- rep(mean(y[pred.cells]), length(y))
      lambda.max <- 0
      lambda.min <- 0
      sd.cv <- 0
    } else {
      mu <- c(glmnet::predict.cv.glmnet(cv, newx = x, s = "lambda.min",
                                        type="response"))
      lambda.max <- cv$lambda[1]
      lambda.min <- cv$lambda.min
      min.ind <- which(cv$lambda == cv$lambda.min)
      sd.cv <- (cv$cvm[1] - cv$cvm[min.ind]) / cv$cvsd[min.ind]
    }
  } else {
    lambda.seq <- c(exp(seq(log(lambda.max), log(lambda.min), by = -0.2)),
                    lambda.min)
    cv <- tryCatch(
      suppressWarnings(glmnet::glmnet(x[pred.cells, ], y[pred.cells],
                                      family="poisson", dfmax = 300,
                                      lambda = lambda.seq)),
      error = function(cond) {
        return(NA)
      }
    )
    if (length(cv) == 1) {
      mu <- rep(mean(y[pred.cells]), length(y))
      lambda.max <- 0
      lambda.min <- 0
      sd.cv <- 0
    } else {
      mu <- exp(c(glmnet::predict.glmnet(cv, newx = x, s = lambda.min,
                                     type="response")))
      sd.cv <- NA
    }
  }

  return(list(mu, lambda.max, lambda.min, sd.cv))
}

