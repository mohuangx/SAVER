#' Combines SAVER
#'
#' Combines SAVER objects
#'
#' If SAVER was applied to a dataset for chunks of genes (by specifying
#' \code{pred.genes} and \code{pred.genes.only = TRUE}), this function combines
#' the individual SAVER objects into one SAVER object.
#'
#' @param saver.list List of SAVER objects
#'
#' @return A combined SAVER object
#'
#' @examples
#' data("linnarsson")
#' 
#' \dontrun{
#' a <- saver(linnarsson, pred.genes = 1:5, pred.genes.only = TRUE)
#' b <- saver(linnarsson, pred.genes = 6:10, pred.genes.only = TRUE)
#' ab <- combine.saver(list(a, b))
#' }
#'
#' @export
#'

combine.saver <- function(saver.list) {
  if (!is.null(saver.list[[1]]$alpha)) {
    return(combine.saver.old(saver.list))
  }
  est <- do.call(rbind, lapply(saver.list, `[[`, 1))
  se <- do.call(rbind, lapply(saver.list, `[[`, 2))
  info <- vector("list", 10)
  names(info) <- c("size.factor", "maxcor", "lambda.max", "lambda.min",
                   "sd.cv", "pred.time", "var.time", "cutoff",
                   "lambda.coefs", "total.time")
  info[[1]] <- saver.list[[1]]$info$size.factor
  info.list <- lapply(saver.list, `[[`, 3)
  for (i in 2:9) {
    info[[i]] <- do.call(c, lapply(info.list, `[[`, i))
  }
  info[[10]] <- do.call(sum, (lapply(info.list, `[[`, 10)))
  out <- list(estimate = est, se = se, info = info)
  class(out) <- "saver"
  out
}
