
chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))
