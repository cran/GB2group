scale.gb2 <- function(theta, m) {
  a <- theta[1]
  p <- theta[2]
  q <- theta[3]
  parb <- m/(beta(p + 1 / a, q - 1 / a) / beta(p, q))
  return(parb)
}

scale.ln <- function(theta, m) {
  s <- theta[1]
  parb <- log(m) - s^2 / 2
  return(parb)
}

moment.gb2 <- function(theta, r) {
  a <- theta[1]
  b <- theta[2]
  p <- theta[3]
  q <- theta[4]
  if(q > r/a){
    m.gb2 <- b^r*beta(p + r / a, q - r / a) / beta(p, q)
  }
  if(q <= r/a){
    m.gb2 <- sum(rgb2(10^4, a, b, p, q)^r)/10^4
  }
  return(m.gb2)
}

moment.ln <- function(theta, r) {
  s <- theta[1]
  mu <- theta[2]
  m.ln <- exp(mu * r + r^2 * s^2 / 2)
  return(m.ln)
}
