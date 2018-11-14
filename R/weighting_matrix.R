#' @importFrom stats qlnorm
weight.mat.gb2 <- function(theta, pr) {
  a <- theta[1]
  b <- theta[2]
  p <- theta[3]
  q <- theta[4]
  if(q > 2/a) {
    e1 <- moment.gb2(theta, 1)
    e2 <- moment.gb2(theta, 2)
    w.cop  <- function(x, y) {
      e2*dmoment.gb2(theta, qgb2(x, a, b, p, q), 2) + (x * qgb2(x, a, b, p, q) - e1 * lc.gb2(theta, x)) *
        (qgb2(y, a, b, p, q) - y * qgb2(y, a, b, p, q) + e1 * lc.gb2(theta, y)) - qgb2(x, a, b, p, q) * e1 * lc.gb2(theta, x)
    }
    w.mat <- outer(c(pr, 0.999999999), c(pr, 0.999999999), w.cop) # Quantile of order p = 0.999999999 is used to approximate the maximum
    w.mat[lower.tri(w.mat)] <- t(w.mat)[lower.tri(w.mat)]
    y.cop <- matrix(cbind(diag(1/e1, length(pr)), -lc.gb2(theta, pr) / e1), length(pr), length(pr) + 1)
    return(solve(y.cop %*% w.mat %*% t(y.cop)))
  }
  else {
    stop('q <= 2/a')
  }
}

weight.mat.da <- function(theta, pr) {
  a <- theta[1]
  b <- theta[2]
  p <- theta[3]
  q <- 1
  theta <- c(a, b, p, q)
  if(q > 2/a) {
    e1 <- moment.gb2(theta, 1)
    e2 <- moment.gb2(theta, 2)
    w.cop  <- function(x, y) {
      e2*dmoment.gb2(theta, qgb2(x, a, b, p, q), 2) + (x * qgb2(x, a, b, p, q) - e1 * lc.gb2(theta, x)) *
        (qgb2(y, a, b, p, q) - y * qgb2(y, a, b, p, q) + e1 * lc.gb2(theta, y)) - qgb2(x, a, b, p, q) * e1 * lc.gb2(theta, x)
    }
    w.mat <- outer(c(pr, 0.999999999), c(pr, 0.999999999), w.cop) # Quantile of order p = 0.999999999 is used to approximate the maximum
    w.mat[lower.tri(w.mat)] <- t(w.mat)[lower.tri(w.mat)]
    y.cop <- matrix(cbind(diag(1/e1, length(pr)), -lc.gb2(theta, pr) / e1), length(pr), length(pr) + 1)
    return(solve(y.cop %*% w.mat %*% t(y.cop)))
  }
  else {
    stop('1 <= 2/a')
  }
}

weight.mat.b2 <- function(theta, pr) {
  a <- 1
  b <- theta[1]
  p <- theta[2]
  q <- theta[3]
  theta <- c(a, b, p, q)
  if(q > 2) {
    e1 <- moment.gb2(theta, 1)
    e2 <- moment.gb2(theta, 2)
    w.cop  <- function(x, y) {
      e2*dmoment.gb2(theta, qgb2(x, a, b, p, q), 2) + (x * qgb2(x, a, b, p, q) - e1 * lc.gb2(theta, x)) *
        (qgb2(y, a, b, p, q) - y * qgb2(y, a, b, p, q) + e1 * lc.gb2(theta, y)) - qgb2(x, a, b, p, q) * e1 * lc.gb2(theta, x)
    }
    w.mat <- outer(c(pr, 0.999999999), c(pr, 0.999999999), w.cop) # Quantile of order p = 0.999999999 is used to approximate the maximum
    w.mat[lower.tri(w.mat)] <- t(w.mat)[lower.tri(w.mat)]
    y.cop <- matrix(cbind(diag(1/e1, length(pr)), -lc.gb2(theta, pr) / e1), length(pr), length(pr) + 1)
    return(solve(y.cop %*% w.mat %*% t(y.cop)))
  }
}

weight.mat.sm <- function(theta, pr) {
  a <- theta[1]
  b <- theta[2]
  p <- 1
  q <- theta[3]
  theta <- c(a, b, p, q)
  if(q > 2/a) {
    e1 <- moment.gb2(theta, 1)
    e2 <- moment.gb2(theta, 2)
    w.cop  <- function(x, y) {
      e2*dmoment.gb2(theta, qgb2(x, a, b, p, q), 2) + (x * qgb2(x, a, b, p, q) - e1 * lc.gb2(theta, x)) *
        (qgb2(y, a, b, p, q) - y * qgb2(y, a, b, p, q) + e1 * lc.gb2(theta, y)) - qgb2(x, a, b, p, q) * e1 * lc.gb2(theta, x)
    }
    w.mat <- outer(c(pr, 0.999999999), c(pr, 0.999999999), w.cop) # Quantile of order p = 0.999999999 is used to approximate the maximum
    w.mat[lower.tri(w.mat)] <- t(w.mat)[lower.tri(w.mat)]
    y.cop <- matrix(cbind(diag(1/e1, length(pr)), -lc.gb2(theta, pr) / e1), length(pr), length(pr) + 1)
    return(solve(y.cop %*% w.mat %*% t(y.cop)))
  }
  else {
    stop('q <= 2/a')
  }
}

weight.mat.f <- function(theta, pr) {
  a <- theta[1]
  b <- theta[2]
  p <- 1
  q <- 1
  theta <- c(a, b, p, q)
  if(1 > 2/a) {
    e1 <- moment.gb2(theta, 1)
    e2 <- moment.gb2(theta, 2)
    w.cop  <- function(x, y) {
      e2*dmoment.gb2(theta, qgb2(x, a, b, p, q), 2) + (x * qgb2(x, a, b, p, q) - e1 * lc.gb2(theta, x)) *
        (qgb2(y, a, b, p, q) - y * qgb2(y, a, b, p, q) + e1 * lc.gb2(theta, y)) - qgb2(x, a, b, p, q) * e1 * lc.gb2(theta, x)
    }
    w.mat <- outer(c(pr, 0.999999999), c(pr, 0.999999999), w.cop) # Quantile of order p = 0.999999999 is used to approximate the maximum
    w.mat[lower.tri(w.mat)] <- t(w.mat)[lower.tri(w.mat)]
    y.cop <- matrix(cbind(diag(1/e1, length(pr)), -lc.gb2(theta, pr) / e1), length(pr), length(pr) + 1)
    return(solve(y.cop %*% w.mat %*% t(y.cop)))
  }
  else {
    stop('1 <= 2/a')
  }
}

weight.mat.ln <- function(theta, pr) {
  s <- theta[1]
  mu <- theta[2]
  e1 <- moment.ln(theta, 1)
  e2 <- moment.ln(theta, 2)
  w.cop  <- function(x, y) {
    e2 * dmoment.ln(theta, qlnorm(x, mu, s), 2) + (x * qlnorm(x, mu, s) - e1 * lc.ln(theta, x)) *
      (qlnorm(y, mu, s) - y * qlnorm(y, mu, s) + e1 * lc.ln(theta, y)) - qlnorm(x, mu, s) * e1 * lc.ln(theta, x)
  }
  w.mat <- outer(c(pr, 0.999999999), c(pr, 0.999999999), w.cop)
  w.mat[lower.tri(w.mat)] <- t(w.mat)[lower.tri(w.mat)]
  y.cop <- matrix(cbind(diag(1 / e1, length(pr)), -lc.ln(theta, pr) / e1), length(pr), length(pr) + 1)
  return(solve(y.cop %*% w.mat %*% t(y.cop)))
}

