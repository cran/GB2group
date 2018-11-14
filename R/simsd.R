#' @importFrom ineq Lc
#' @importFrom stats sd rlnorm

simsd.gb2 <- function(x, theta, N, nrep, se.scale = F) {
  if(is.na(theta[2])) {
    a <- theta[1]
    b <- 1
    p <- theta[3]
    q <- theta[4]
  }
  else {
    a <- theta[1]
    b <- theta[2]
    p <- theta[3]
    q <- theta[4]
  }
  cprob <- as.vector(x[!is.na(x)])
  cprob <- as.numeric(cprob)/sum(cprob)
  cprob <- cumsum(cprob)[-length(cprob)]
  par.sim <- matrix(NA, nrep, 4)
  gmm.sc <- rep(NA, nrep)
  for(i in 1:nrep) {
    sim.sam <- rgb2(N, a, b, p, q)
    share <- Lc(sim.sam)$L[cprob * N + 1]
    regress <- try(suppressWarnings(nlsLM(share ~ (lc.gb2(c(A, 1, P, Q), cprob)), algorithm = "port",
      start = list(A = a, P = p, Q = q), lower = c(0, 0, 0), control = nls.lm.control(maxiter = 1000))),
      silent = TRUE)
    if('try-error'%in%class(regress)) next
    else {
      par.sim[i, c(1, 3, 4)] <- coef(regress)
      if(is.na(theta[2])) next
      else{
        sim.mean <- mean(sim.sam)
        par.sim[i, 2] <- scale.gb2(coef(regress), sim.mean)
        if(se.scale == TRUE ){
          regress <- try(suppressWarnings(opt.gmm.gb2(cprob, share, init.est = par.sim[i, ], cons.est = par.sim[i, ])))
          if(!'try-error'%in%class(regress$opt1)) {
            gmm.sc[i] <- scale.gb2(regress$opt1$par[-2], sim.mean)
          }
        }
      }
    }
  }
  sd.scale <- NA
  if(se.scale == TRUE) {
    sd.scale <- sd(gmm.sc, na.rm = T)
  }
  nls.se <- c(sd(par.sim[, 1], na.rm = T), sd(par.sim[, 2], na.rm = T), sd(par.sim[, 3], na.rm = T), sd(par.sim[, 4], na.rm = T))
  return(list(nls.se = nls.se, sd.scale = sd.scale))
}

simsd.da <- function(x, theta, N, nrep, se.scale = FALSE) {
  if(is.na(theta[2])) {
    a <- theta[1]
    b <- 1
    p <- theta[3]
    q <- 1
  }
  else {
    a <- theta[1]
    b <- theta[2]
    p <- theta[3]
    q <- 1
  }
  cprob <- as.vector(x[!is.na(x)])
  cprob <- as.numeric(cprob)/sum(cprob)
  cprob <- cumsum(cprob)[-length(cprob)]
  par.sim <- matrix(NA, nrep, 3)
  gmm.sc <- rep(NA, nrep)
  for(i in 1:nrep) {
    sim.sam <- rgb2(N, a, b, p, q)
    share <- Lc(sim.sam)$L[cprob * N + 1]
    regress <- try(suppressWarnings(nlsLM(share ~ (lc.gb2(c(A, 1, P, 1), cprob)), algorithm = "port",
      start = list(A = a, P = p), lower = c(0, 0), control = nls.lm.control(maxiter = 1000))),
      silent = TRUE)
    if('try-error'%in%class(regress)) next
    else {
      par.sim[i, c(1, 3)] <- coef(regress)
      if(is.na(theta[2])) next
      else{
        sim.mean <- mean(sim.sam)
        par.sim[i, 2] <- scale.gb2(c(coef(regress), 1), sim.mean)
        if(se.scale == TRUE & (1 > 2 / par.sim[i, 1])){
          regress <- try(opt.gmm.da(cprob, share, init.est = par.sim[i, ], cons.est = par.sim[i, ]))
          if('try-error'%in%class(regress$opt1)) next
          else{
            gmm.sc[i] <- scale.gb2(c(regress$opt1$par[-2], 1), sim.mean)
          }
        }
      }
    }
  }
  sd.scale <- NA
  if(se.scale == TRUE) {
    sd.scale <- sd(gmm.sc, na.rm = T)
  }
  nls.se <- c(sd(par.sim[, 1], na.rm = T), sd(par.sim[, 2], na.rm = T), sd(par.sim[, 3], na.rm = T))
  return(list(nls.se = nls.se, sd.scale = sd.scale))
}

simsd.b2 <- function(x, theta, N, nrep, se.scale = FALSE) {
  if(is.na(theta[1])) {
    a <- 1
    b <- 1
    p <- theta[2]
    q <- theta[3]
  }
  else {
    a <- 1
    b <- theta[1]
    p <- theta[2]
    q <- theta[3]
  }
  cprob <- as.vector(x[!is.na(x)])
  cprob <- as.numeric(cprob)/sum(cprob)
  cprob <- cumsum(cprob)[-length(cprob)]
  par.sim <- matrix(NA, nrep, 3)
  gmm.sc <- rep(NA, nrep)
  for(i in 1:nrep) {
    sim.sam <- rgb2(N, a, b, p, q)
    share <- Lc(sim.sam)$L[cprob * N + 1]
    regress <- try(suppressWarnings(nlsLM(share ~ (lc.gb2(c(1, 1, P, Q), cprob)), algorithm = "port",
      start = list(P = p, Q = q), lower = c(0, 0), control = nls.lm.control(maxiter = 1000))),
      silent = TRUE)
    if('try-error'%in%class(regress)) next
    else {
      par.sim[i, c(2, 3)] <- coef(regress)
      if(is.na(theta[1])) next
      else{
        sim.mean <- mean(sim.sam)
        par.sim[i, 1] <- scale.gb2(c(1, coef(regress)), sim.mean)
        if(se.scale == TRUE & (par.sim[i, 3] > 2)){
          regress <- try(opt.gmm.b2(cprob, share, init.est = par.sim[i, ], cons.est = par.sim[i, ]))
          if('try-error'%in%class(regress$opt1)) next
          else{
            gmm.sc[i] <- scale.gb2(c(1, regress$opt1$par[-1]), sim.mean)
          }
        }
      }
    }
  }
  sd.scale <- NA
  if(se.scale == TRUE) {
    sd.scale <- sd(gmm.sc, na.rm = T)
  }
  nls.se <- c(sd(par.sim[, 1], na.rm = T), sd(par.sim[, 2], na.rm = T), sd(par.sim[, 3], na.rm = T))
  return(list(nls.se = nls.se, sd.scale = sd.scale))
}


simsd.sm <- function(x, theta, N, nrep, se.scale = FALSE) {
  if(is.na(theta[2])) {
    a <- theta[1]
    b <- 1
    p <- 1
    q <- theta[3]
  }
  else {
    a <- theta[1]
    b <- theta[2]
    p <- 1
    q <- theta[3]
  }
  cprob <- as.vector(x[!is.na(x)])
  cprob <- as.numeric(cprob)/sum(cprob)
  cprob <- cumsum(cprob)[-length(cprob)]
  par.sim <- matrix(NA, nrep, 3)
  gmm.sc <- rep(NA, nrep)
  for(i in 1:nrep) {
    sim.sam <- rgb2(N, a, b, p, q)
    share <- Lc(sim.sam)$L[cprob * N + 1]
    regress <- try(suppressWarnings(nlsLM(share ~ (lc.gb2(c(A, 1, 1, Q), cprob)), algorithm = "port",
      start = list(A = a, Q = q), lower = c(0, 0), control = nls.lm.control(maxiter = 1000))),
      silent = TRUE)
    if('try-error'%in%class(regress)) next
    else {
      par.sim[i, c(1, 3)] <- coef(regress)
      if(is.na(theta[2])) next
      else{
        sim.mean <- mean(sim.sam)
        par.sim[i, 2] <- scale.gb2(c(coef(regress)[1], 1, coef(regress)[2]), sim.mean)
        if(se.scale == TRUE & (par.sim[i, 3] > 2 / par.sim[i, 1])){
          regress <- try(opt.gmm.sm(cprob, share, init.est = par.sim[i, ], cons.est = par.sim[i, ]))
          if('try-error'%in%class(regress$opt1)) next
          else{
            gmm.sc[i] <- scale.gb2(c(regress$opt1$par[1], 1, regress$opt1$par[3]), sim.mean)
          }
        }
      }
    }
  }
  sd.scale <- NA
  if(se.scale == TRUE) {
    sd.scale <- sd(gmm.sc, na.rm = T)
  }
  nls.se <- c(sd(par.sim[, 1], na.rm = T), sd(par.sim[, 2], na.rm = T), sd(par.sim[, 3], na.rm = T))
  return(list(nls.se = nls.se, sd.scale = sd.scale))
}

simsd.f <- function(x, theta, N, nrep, se.scale = FALSE) {
  if(is.na(theta[2])) {
    a <- theta[1]
    b <- 1
    p <- 1
    q <- 1
  }
  else {
    a <- theta[1]
    b <- theta[2]
    p <- 1
    q <- 1
  }
  cprob <- as.vector(x[!is.na(x)])
  cprob <- as.numeric(cprob)/sum(cprob)
  cprob <- cumsum(cprob)[-length(cprob)]
  par.sim <- matrix(NA, nrep, 2)
  gmm.sc <- rep(NA, nrep)
  for(i in 1:nrep) {
    sim.sam <- rgb2(N, a, b, p, q)
    share <- Lc(sim.sam)$L[cprob * N + 1]
    regress <- try(suppressWarnings(nlsLM(share ~ (lc.gb2(c(A, 1, 1, 1), cprob)), algorithm = "port",
      start = list(A = a), lower = 0, control = nls.lm.control(maxiter = 1000))),
      silent = TRUE)
    if('try-error'%in%class(regress)) next
    else {
      par.sim[i, 1] <- coef(regress)
      if(is.na(theta[2])) next
      else{
        sim.mean <- mean(sim.sam)
        par.sim[i, 2] <- scale.gb2(c(coef(regress)[1], 1, 1), sim.mean)
        if(se.scale == TRUE & (1 > 2 / par.sim[i, 1])){
          regress <- try(opt.gmm.f(cprob, share, init.est = par.sim[i, ], cons.est = par.sim[i, ]))
          if('try-error'%in%class(regress$opt1)) next
          else{
            gmm.sc[i] <- scale.gb2(c(regress$opt1$par[1], 1, 1), sim.mean)
          }
        }
      }
    }
  }
  sd.scale <- NA
  if(se.scale == TRUE) {
    sd.scale <- sd(gmm.sc, na.rm = T)
  }
  nls.se <- c(sd(par.sim[, 1], na.rm = T), sd(par.sim[, 2], na.rm = T))
  return(list(nls.se = nls.se, sd.scale = sd.scale))
}

simsd.ln <- function(x, theta, N, nrep, se.scale = FALSE) {
  if(is.na(theta[2])) {
    s <- theta[1]
    mu <- 1
  }
  else {
    s <- theta[1]
    mu <- theta[2]
  }
  cprob <- as.vector(x[!is.na(x)])
  cprob <- as.numeric(cprob)/sum(cprob)
  cprob <- cumsum(cprob)[-length(cprob)]
  par.sim <- matrix(NA, nrep, 2)
  gmm.sc <- rep(NA, nrep)
  for(i in 1:nrep) {
    sim.sam <- rlnorm(N, mu, s)
    share <- Lc(sim.sam)$L[cprob * N + 1]
    regress <- try(suppressWarnings(nlsLM(share ~ (lc.ln(S, cprob)), algorithm ="port", start = list(S = s), lower = 0, control = nls.lm.control(maxiter=1000))),
      silent = TRUE)
    if('try-error'%in%class(regress)) next
    else {
      par.sim[i, 1] <- coef(regress)
      if(is.na(theta[2])) next
      else{
        sim.mean <- mean(sim.sam)
        par.sim[i, 2] <- scale.ln(coef(regress)[1], sim.mean)
        if(se.scale == TRUE){
          regress <- try(opt.gmm.ln(cprob, share, init.est = par.sim[i, ], cons.est = par.sim[i, ]))
          if('try-error'%in%class(regress$opt1)) next
          else{
            gmm.sc[i] <- scale.ln(regress$opt1$par[1], sim.mean)
          }
        }
      }
    }
  }
  sd.scale <- NA
  if(se.scale == TRUE) {
    sd.scale <- sd(gmm.sc, na.rm = T)
  }
  nls.se <- c(sd(par.sim[, 1], na.rm = T), sd(par.sim[, 2], na.rm = T))
  return(list(nls.se = nls.se, sd.scale = sd.scale))
}