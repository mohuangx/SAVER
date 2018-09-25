#' Output prior predictions
#' 
#' Outputs prior predictions
#' 
#' This function outputs prior mean predictions \eqn{\mu} used in fitting the 
#' SAVER model.
#' 
#' @param x Original count matrix.
#' 
#' @param saver.obj SAVER output.
#' 
#' @return A matrix of prior mean predictions
#' 
#' @examples 
#' data("linnarsson")
#' 
#' saver1 <- saver(linnarsson, pred.genes = 1:5, pred.genes.only = TRUE)
#' 
#' mu <- get.mu(linnarsson[1:5, ], saver1)
#' mu[1:5, 1:5]
#' 

get.mu <- function(x, saver.obj) {
  mu <- (saver.obj$estimate^2/saver.obj$se^2 - x)/
    sweep(saver.obj$estimate/saver.obj$se^2, 2, saver.obj$info$size.factor, "-")
  round(mu, 3)
}