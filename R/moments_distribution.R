dmoment.gb2 <- function(theta, x, r) {
  a <- theta[1]
  b <- theta[2]
  p <- theta[3]
  q <- theta[4]
  dm.gb2 <- pgb2(x, a, b, p + r/a, q - r/a)
  return(dm.gb2)
}

dmoment.ln <- function(theta, x, r) {
  s <- theta[1]
  mu <- theta[2]
  dm.gb2 <- pnorm((log(x) - mu - r * s^2)/s)
  return(dm.gb2)
}

