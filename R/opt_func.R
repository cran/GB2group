opt.gmm.gb2 <- function (x, y, init.est, cons.est, est.method = 1) {
  minfun <- function(theta, x, y) {
    moments.gmm <- matrix((lc.gb2(theta, x) - y), length(x), 1)
    min.gmm <- t(moments.gmm) %*% weight.mat.gb2(cons.est, x) %*% moments.gmm
    return(min.gmm)
  }
  opt1 <- try(optim(init.est, minfun, gr = NULL, x, y, method = "BFGS", control = list(parscale = init.est, pgtol = 1e-08)))
  if (est.method != 2) {
    return(list(opt1 = opt1))
  }
  if (est.method == 2 |'try-error'%in%class(opt1)) {
    opt2 <- optim(init.est, minfun, gr = NULL, x, y, method = "L-BFGS-B", lower = 0, control = list(parscale = init.est, pgtol = 1e-08))
    return(list(opt2 = opt2))
  }
}

opt.gmm.da <- function (x, y, init.est, cons.est, est.method = 1) {
  minfun <- function(theta, x, y) {
    theta <- c(theta, 1)
    moments.gmm <- matrix((lc.gb2(theta, x) - y), length(x), 1)
    min.gmm <- t(moments.gmm) %*% weight.mat.da(cons.est, x) %*% moments.gmm
    return(min.gmm)
  }
  opt1 <- try(optim(init.est, minfun, gr = NULL, x, y, method = "BFGS", control = list(parscale = init.est, pgtol = 1e-08)))
  if (est.method != 2) {
    return(list(opt1 = opt1))
  }
  if (est.method == 2 |'try-error'%in%class(opt1)) {
    opt2 <- optim(init.est, minfun, gr = NULL, x, y, method = "L-BFGS-B", lower = 0, control = list(parscale = init.est, pgtol = 1e-08))
    return(list(opt2 = opt2))
  }
}

opt.gmm.b2 <- function (x, y, init.est, cons.est, est.method = 1) {
  minfun <- function(theta, x, y) {
    theta <- c(1, theta)
    moments.gmm <- matrix((lc.gb2(theta, x) - y), length(x), 1)
    min.gmm <- t(moments.gmm) %*% weight.mat.b2(cons.est, x) %*% moments.gmm
    return(min.gmm)
  }
  opt1 <- try(optim(init.est, minfun, gr = NULL, x, y, method = "BFGS", control = list(parscale = init.est, pgtol = 1e-08)))
  if (est.method != 2) {
    return(list(opt1 = opt1))
  }
  if (est.method == 2 |'try-error'%in%class(opt1)) {
    opt2 <- optim(init.est, minfun, gr = NULL, x, y, method = "L-BFGS-B", lower = 0, control = list(parscale = init.est, pgtol = 1e-08))
    return(list(opt2 = opt2))
  }
}

opt.gmm.sm <- function (x, y, init.est, cons.est, est.method = 1) {
  minfun <- function(theta, x, y) {
    theta <- c(theta[c(1, 2)], 1, theta[3])
    moments.gmm <- matrix((lc.gb2(theta, x) - y), length(x), 1)
    min.gmm <- t(moments.gmm) %*% weight.mat.sm(cons.est, x) %*% moments.gmm
    return(min.gmm)
  }
  opt1 <- try(optim(init.est, minfun, gr = NULL, x, y, method = "BFGS", control = list(parscale = init.est, pgtol = 1e-08)))
  if (est.method != 2) {
    return(list(opt1 = opt1))
  }
  if (est.method == 2 |'try-error'%in%class(opt1)) {
    opt2 <- optim(init.est, minfun, gr = NULL, x, y, method = "L-BFGS-B", lower = 0, control = list(parscale = init.est, pgtol = 1e-08))
    return(list(opt2 = opt2))
  }
}

opt.gmm.f <- function (x, y, init.est, cons.est, est.method = 1) {
  minfun <- function(theta, x, y) {
    theta <- c(theta, 1, 1)
    moments.gmm <- matrix((lc.gb2(theta, x) - y), length(x), 1)
    min.gmm <- t(moments.gmm) %*% weight.mat.f(cons.est, x) %*% moments.gmm
    return(min.gmm)
  }
  opt1 <- try(optim(init.est, minfun, gr = NULL, x, y, method = "BFGS", control = list(parscale = init.est, pgtol = 1e-08)))
  if (est.method != 2) {
    return(list(opt1 = opt1))
  }
  if (est.method == 2 |'try-error'%in%class(opt1)) {
    opt2 <- optim(init.est, minfun, gr = NULL, x, y, method = "L-BFGS-B", lower = 0, control = list(parscale = init.est, pgtol = 1e-08))
    return(list(opt2 = opt2))
  }
}

opt.gmm.ln <- function (x, y, init.est, cons.est, est.method = 1, hess = FALSE) {
  minfun <- function(theta, x, y) {
    moments.gmm <- matrix((lc.ln(theta, x) - y), length(x), 1)
    min.gmm <- t(moments.gmm) %*% weight.mat.ln(cons.est, x) %*% moments.gmm
    return(min.gmm)
  }
  opt1 <- optim(init.est, minfun, gr = NULL, x, y, method = "BFGS", control = list(parscale = init.est, pgtol = 1e-08), hessian = hess)
  if (est.method != 2)
    return(list(opt1 = opt1))
  if (est.method == 2 |'try-error'%in%class(opt1)) {
    opt2 <- optim(init.est, minfun, gr = NULL, x, y, method = "L-BFGS-B", lower = 0, control = list(parscale = init.est, pgtol = 1e-08), hessian = hess)
    return(list(opt2 = opt2))
  }
}

