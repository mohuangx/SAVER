
expr.predict <- function(x, y, dfmax = 300, nfolds = 5, seed = NULL) {
  if (!is.null(seed))
    set.seed(seed)
  cv <- tryCatch(
    suppressWarnings(glmnet::cv.glmnet(x, y, family="poisson", dfmax = dfmax,
                               nfolds = nfolds)),
    error = function(cond) {
      return(NA)
    }
  )
  if (length(cv) == 1) {
    mu <- rep(mean(y), length(y))
  } else {
    mu <- c(glmnet::predict.cv.glmnet(cv, newx = x, s = "lambda.min",
                                      type="response"))
  }
  return(mu)
}
