
gr.gb2 <- function(theta) {
  a <- theta[1]
  b <- theta[2]
  p <- theta[3]
  q <- theta[4]
  pr <- theta[5]
  t <- qgb2(pr, a, b, p, q)
  lcgb2 <- pgb2(t, a, b, p + 1 / a, q - 1 / a)
  return(-lcgb2)
}

gr.da <- function(theta) {
  a <- theta[1]
  b <- theta[2]
  p <- theta[3]
  q <- 1
  pr <- theta[4]
  t <- qgb2(pr, a, b, p, q)
  lcgb2 <- pgb2(t, a, b, p + 1 / a, q - 1 / a)
  return(-lcgb2)
}

gr.b2 <- function(theta) {
  a <- 1
  b <- theta[1]
  p <- theta[2]
  q <- theta[3]
  pr <- theta[4]
  t <- qgb2(pr, a, b, p, q)
  lcgb2 <- pgb2(t, a, b, p + 1 / a, q - 1 / a)
  return(-lcgb2)
}

gr.sm <- function(theta) {
  a <- theta[1]
  b <- theta[2]
  p <- 1
  q <- theta[3]
  pr <- theta[4]
  t <- qgb2(pr, a, b, p, q)
  lcgb2 <- pgb2(t, a, b, p + 1 / a, q - 1 / a)
  return(-lcgb2)
}

gr.f <- function(theta) {
  a <- theta[1]
  b <- theta[2]
  p <- 1
  q <- 1
  pr <- theta[3]
  t <- qgb2(pr, a, b, p, q)
  lcgb2 <- pgb2(t, a, b, p + 1 / a, q - 1 / a)
  return(-lcgb2)
}

gr.ln <- function(theta) {
  s <- theta[1]
  mu <- theta[2]
  pr <- theta[3]
  lcln <- pnorm(qnorm(pr) - s)
  return(-lcln)
}

