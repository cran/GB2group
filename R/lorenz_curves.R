#' @importFrom GB2 qgb2 pgb2
#' @importFrom stats pnorm qnorm
lc.gb2 <- function(theta, pr) {
  a <- theta[1]
  b <- theta[2]
  p <- theta[3]
  q <- theta[4]
  t <- qgb2(pr, a, b, p, q)
  lcgb2 <- pgb2(t, a, b, p + 1 / a, q - 1 / a)
  return(lcgb2)
}

lc.ln <- function(theta, pr) {
  s <- theta[1]
  mu <- theta[2]
  lcln <- pnorm(qnorm(pr) - s)
  return(lcln)
}
