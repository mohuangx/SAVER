
expr.predict <- function(x, y, dfmax = 300, nfolds = 5) {
  cv <- tryCatch(
    suppressWarnings(cv.glmnet(x, y, family="poisson", dfmax = dfmax,
                               nfolds = nfolds)),
    error = function(cond) {
      return(NA)
    }
  )
  if (length(cv) == 1) {
    mu <- rep(mean(y), length(y))
  } else { # return prediction
    mu <- c(predict(cv, newx = x, s = "lambda.min", type="response"))
  }
  return(mu)
}
