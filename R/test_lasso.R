library(microbenchmark)


microbenchmark(colSums(x), apply(x, 2, sum), .colSums(x, nrow(x), ncol(x)))
microbenchmark(sweep(x[1:500, ], 2, sf, "/"))

library(penalized)
library(glmnet)



g1 <- function(i) {
  glmnet(x.est[, -i], x[i, ]/sf, family = "poisson", dfmax = 300,
         lambda = c(exp(seq(log(lambda.max), log(lambda.min), by = -0.1)),
                    lambda.min))
}

g2 <- function(i) {
  glmnet(x.est[, -i], x[i, ]/sf, family = "poisson", dfmax = 300,
         lambda = lambda.min)
}

g3 <- function(i) {
  glmnet(x.est[, -i], x[i, ]/sf, family = "poisson", dfmax = 300,
         lambda = c(exp(seq(log(lambda.max), log(lambda.min), by = -0.2)),
                    lambda.min))
}

g4 <- function(i) {
  glmnet(x.est[, -i], x[i, ]/sf, family = "poisson", dfmax = 300,
         lambda = c(exp(seq(log(lambda.max), log(lambda.min), by = -0.5)),
                    lambda.min))
}



i <- 11

r <- calc.maxcor(x.est[, -i], x[i, ]/sf)

lambda.max <- r*sd(x[i, ]/sf)*sqrt(499/500)

lambda.min <- est.lambda.min(x[i, ], r, sf, 
                             fit$coefficients)

cv1 <- cv.glmnet(x.est[, -i], x[i, ]/sf, family = "poisson", dfmax = 300)

plot(cv1)
abline(v = log(lambda.min))

microbenchmark(g1(i), g2(i), g3(i), g4(i), times = 20)
