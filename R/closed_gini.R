gini.b2 <- function(theta) {
  p <- theta[1]
  q <- theta[2]
  2 * beta(2 * p, 2 * q - 1) / (p * beta(p, q)^2)
}

gini.sm <- function(theta) {
  a <- theta[1]
  q <- theta[2]
  1 - gamma(q) * gamma(2 * q - 1 / a) / (gamma(q - 1 / a) * gamma(2 * q))
}

gini.d <- function(theta) {
  a <- theta[1]
  p <- theta[2]
  gamma(p) * gamma(2 * p + 1 / a) / (gamma(p + 1 / a)*gamma(2 * p)) - 1
}

gini.ln <- function(theta) {
  s <- theta[1]
  2 * pnorm(s / 2^0.5) - 1
}
