#' @importFrom GB2 rgb2 gini.gb2 gb2.gini
#' @importFrom ineq Gini
simgini.gb2 <-  function(theta, size = 10^6) {
  if(length(theta) == 3 & sum(is.na(theta)) == 0) {
    a <- theta[1]
    b <- 1
    p <- theta[2]
    q <- theta[3]
  }
  if((length(theta) == 4 & is.na(theta[2])) & sum(is.na(theta)) == 1) {
    a <- theta[1]
    b <- 1
    p <- theta[3]
    q <- theta[4]
  }
  if(length(theta) == 4 & sum(is.na(theta)) == 0) {
    a <- theta[1]
    b <- theta[2]
    p <- theta[3]
    q <- theta[4]
  }
  gini.t <- gini.gb2(a, p, q)
  if (is.nan(gini.t) | is.na(gini.t) | gini.t == 1) {
    sim <- rgb2(size, a, b, p, q)
    sim <- sim[is.finite(sim)]
    gini.t <- Gini(sim, na.rm = TRUE)
  }
  return(gini.t)
}

simgini.da <-  function(theta, size = 10^6) {
  if(length(theta) == 2 & sum(is.na(theta)) == 0) {
    a <- theta[1]
    b <- 1
    p <- theta[2]
    q <- 1
  }
  if((length(theta) == 3 & is.na(theta[2])) & sum(is.na(theta)) == 1) {
    a <- theta[1]
    b <- 1
    p <- theta[3]
    q <- 1
  }
  if(length(theta) == 3 & sum(is.na(theta)) == 0) {
    a <- theta[1]
    b <- theta[2]
    p <- theta[3]
    q <- 1
  }
  gini.t <- gini.d(c(a,p))
  if (is.nan(gini.t)) {
    sim <- rgb2(size, a, b, p, q)
    sim <- sim[is.finite(sim)]
    gini.t <- Gini(sim, na.rm = TRUE)
  }
  return(gini.t)
}

simgini.b2 <-  function(theta, size = 10^6) {
  if(length(theta) == 2 & sum(is.na(theta)) == 0) {
    a <- 1
    b <- 1
    p <- theta[2]
    q <- theta[3]
  }
  if((length(theta) == 3 & is.na(theta[1])) & sum(is.na(theta)) == 1) {
    a <- 1
    b <- 1
    p <- theta[2]
    q <- theta[3]
  }
  if(length(theta) == 3 & sum(is.na(theta)) == 0) {
    a <- 1
    b <- theta[1]
    p <- theta[2]
    q <- theta[3]
  }
  gini.t <- gini.b2(c(p, q))
  if(is.nan(gini.t)){
    sim <- rgb2(size, a, b, p, q)
    sim <- sim[is.finite(sim)]
    gini.t <- Gini(sim, na.rm = TRUE)
  }
  return(gini.t)
}


simgini.sm <-  function(theta, size = 10^6) {
  if(length(theta) == 2 & sum(is.na(theta)) == 0) {
    a <- theta[1]
    b <- 1
    p <- 1
    q <- theta[2]
  }
  if((length(theta) == 3 & is.na(theta[2])) & sum(is.na(theta)) == 1) {
    a <- theta[1]
    b <- 1
    p <- 1
    q <- theta[3]
  }
  if(length(theta) == 3 & sum(is.na(theta)) == 0) {
    a <- theta[1]
    b <- theta[2]
    p <- 1
    q <- theta[3]
  }
  gini.t <- gini.sm(c(a, q))
  if (is.nan(gini.t)) {
    sim <- rgb2(size, a, b, p, q)
    sim <- sim[is.finite(sim)]
    gini.t <- Gini(sim, na.rm = TRUE)
  }
  return(gini.t)
}


simgini.f <-  function(theta, size = 10^6) {
  a <- theta[1]
  gini.t <- 1 / a
  if(gini.t > 1) {
    sim <- rgb2(size, theta[1], 1, 1, 1)
    sim <- sim[is.finite(sim)]
    gini.t <- Gini(sim, na.rm = TRUE)
  }
  return(gini.t)
}

