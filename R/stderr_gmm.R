#' @importFrom numDeriv grad
gmmse.gb2 <- function(theta, pr, N) {
  if(length(theta) != 4 | (length(theta) == 4 & sum(is.na(theta)) != 0) ) {
    stop("Incorrect number of parameters. The GB2 distribution characterized by one scale parameter and three shape paramters")
  }
  derivm <- matrix(NA, length(theta), length(pr))
  for(i in 1:length(pr)){
    gr.lc <- grad(gr.gb2, c(theta, pr[i]), method = "Richardson")
    derivm[, i] <- as.matrix(gr.lc[-length(gr.lc)])
  }
  m1 <- derivm %*% weight.mat.gb2(theta, pr)%*% t(derivm)
  m1 <- m1[-2, ]
  m1 <- m1[, -2]
  mat1 <- solve(m1)
  result <- (diag(mat1)/N)^0.5
  return(c(result[1], NA, result[-1]))
}

gmmse.da <- function(theta, pr, N) {
  if(length(theta) != 3 | (length(theta) == 3 & sum(is.na(theta)) != 0) ) {
    stop("Incorrect number of parameters. The Dagum distribution characterized by one scale parameter and two shape paramters")
  }
  derivm <- matrix(NA, length(theta), length(pr))
  for(i in 1:length(pr)){
    gr.lc <- grad(gr.da, c(theta, pr[i]), method = "Richardson")
    derivm[, i] <- as.matrix(gr.lc[-length(gr.lc)])
  }
  m1 <- derivm %*% weight.mat.da(theta, pr)%*% t(derivm)
  m1 <- m1[-2, ]
  m1 <- m1[, -2]
  mat1 <- solve(m1)
  result <- (diag(mat1)/N)^0.5
  return(c(result[1], NA, result[-1]))
}

gmmse.b2 <- function(theta, pr, N) {
  if(length(theta) != 3 | (length(theta) == 3 & sum(is.na(theta)) != 0) ) {
    stop("Incorrect number of parameters. The B2 distribution characterized by one scale parameter and two shape paramters")
  }
  derivm <- matrix(NA, length(theta), length(pr))
  for(i in 1:length(pr)){
    gr.lc <- grad(gr.b2, c(theta, pr[i]), method = "Richardson")
    derivm[, i] <- as.matrix(gr.lc[-length(gr.lc)])
  }
  m1 <- derivm %*% weight.mat.b2(theta, pr)%*% t(derivm)
  m1 <- m1[-1, ]
  m1 <- m1[, -1]
  mat1 <- solve(m1)
  result <- (diag(mat1)/N)^0.5
  return(c(NA, result))
}

gmmse.sm <- function(theta, pr, N) {
  if(length(theta) != 3 | (length(theta) == 3 & sum(is.na(theta)) != 0) ) {
    stop("Incorrect number of parameters. The Singh-Maddala distribution characterized by one scale parameter and three shape paramters")
  }
  derivm <- matrix(NA, length(theta), length(pr))
  for(i in 1:length(pr)){
    gr.lc <- grad(gr.sm, c(theta, pr[i]), method = "Richardson")
    derivm[, i] <- as.matrix(gr.lc[-length(gr.lc)])
  }
  m1 <- derivm %*% weight.mat.sm(theta, pr)%*% t(derivm)
  m1 <- m1[-2, ]
  m1 <- m1[, -2]
  mat1 <- solve(m1)
  result <- (diag(mat1)/N)^0.5
  return(c(result[1], NA, result[-1]))
}

gmmse.f <- function(theta, pr, N) {
  if(length(theta) != 2 | (length(theta) == 2 & sum(is.na(theta)) != 0) ) {
    stop("Incorrect number of parameters. The Fisk distribution is characterized by one scale parameter and one shape paramter")
  }
  derivm <- matrix(NA, length(theta), length(pr))
  for(i in 1:length(pr)){
    gr.lc <- grad(gr.f, c(theta, pr[i]), method = "Richardson")
    derivm[, i] <- as.matrix(gr.lc[-length(gr.lc)])
  }
  m1 <- derivm %*% weight.mat.f(theta, pr)%*% t(derivm)
  m1 <- m1[-2, ]
  m1 <- m1[-2]
  mat1 <- solve(m1)
  result <- (diag(mat1)/N)^0.5
  return(c(result[1], NA))
}

gmmse.ln <- function(theta, pr, N) {
  if(length(theta) != 2 | (length(theta) == 2 & sum(is.na(theta)) != 0) ) {
    stop("Incorrect number of parameters. The log-normal distribution is characterized by one scale parameter and one shape paramter")
  }
  derivm <- matrix(NA, length(theta), length(pr))
  for(i in 1:length(pr)){
    gr.lc <- grad(gr.ln, c(theta, pr[i]), method = "Richardson")
    derivm[, i] <- as.matrix(gr.lc[-length(gr.lc)])
  }
  m1 <- derivm %*% weight.mat.ln(theta, pr)%*% t(derivm)
  m1 <- m1[-2, ]
  m1 <- m1[-2]
  mat1 <- solve(m1)
  result <- (diag(mat1)/N)^0.5
  return(c(result[1], NA))
}

