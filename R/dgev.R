# Copied from package revdbayes
dgev <- function (x, loc = 0, scale = 1, shape = 0, log = FALSE, m = 1) {
  if (any(scale <= 0)) {
    stop("invalid scale: scale must be positive.")
  }
  max_len <- max(length(x), length(loc), length(scale), length(shape),
                 length(m))
  x <- rep_len(x, max_len)
  loc <- rep_len(loc, max_len)
  scale <- rep_len(scale, max_len)
  shape <- rep_len(shape, max_len)
  m <- rep_len(m, max_len)
  x <- (x - loc) / scale
  xx <- 1 + shape * x
  d <- ifelse(xx < 0, 0,
              ifelse(xx == 0 & shape == -1, 1,
                     ifelse(xx == 0 & shape < -1, Inf,
                            ifelse(abs(shape) > 1e-6,
                                   xx ^ (-(1 + 1 / shape)) * exp(-m * xx ^ (-1/ shape)),
                                   exp(-x + shape * x * (x - 2) / 2 - m * exp(-x + shape * x ^ 2 / 2))))))
  d <- d * m / scale
  if (log) {
    d <- log(d)
  }
  return(d)
}
