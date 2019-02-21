#' Plot of the estimated Lorenz curve and the empirical income shares
#'
#' The function \code{fit.plot} plots the parametric Lorenz curve and the observed income shares used for the estimation of the income distributions belonging to the GB2 family.
#'
#' @param fit  A character string "name" naming the object that contains the estimation of the parametric model for which the Lorenz curve is plotted.
#' @param fit.type specifies the method used to estimate the parametric model. By default, \code{fit.type = 1}, which represents the Lorenz curve estimated by NLS. If \code{fit.type = 2}, the Lorenz curve belongs to the GMM estimation.
#' @param fit.legend If \code{TRUE}, the graph includes a legend indicating the model for which the Lorenz curve is plotted.
#' @param  l.size determines the size of the legend.
#'
#' @return the function \code{fit.plot} returns a graph with the theoretical Lorenz curves of the Generalised Beta of the Second Kind (GB2) family of income distributions and the income shares used for the estimation of these models.
#'
#' @details The function \code{fit.plot} represents the parametric Lorenz curves of some models of the GB2 family. Closed expressions
#'  of the Lorenz curves of these models are provided by Jorda et al. (2018). The parametric model must be estimated
#' before representing the theoretical Lorenz curve. To do so, create an object containing the result of the following functions:
#'  \code{\link{fitgroup.gb2}}, \code{\link{fitgroup.b2}}, \code{\link{fitgroup.da}},
#' \code{\link{fitgroup.sm}}, \code{\link{fitgroup.ln}} or \code{\link{fitgroup.f}}. The name of this object is used (with quotations marks)
#' as the first argument of \code{fit.plot} (see examples below). This function returns a plot with the Lorenz curve of the model
#' estimated by NLS or GMM. More than one fit can be plotted, even when different sets of data are used. The legend indicates the
#' distribution for which the Lorenz curve is represented.
#'
#' @references
#'  Jorda, V., Sarabia, J.M., & Jäntti, M. (2018). Estimation of income inequality from grouped data.
#'  arXiv preprint arXiv:1808.09831.
#'
#' @export
#' @examples
#' fit.ln <- fitgroup.ln(y = c(9, 13, 17, 22, 39), gini.e = 0.29)
#' fit.b2 <- fitgroup.b2(y = c(9, 13, 17, 22, 39), gini.e = 0.29)
#' fit.plot(c("fit.ln", "fit.b2"), fit.legend = TRUE, l.size = 0.8)
#' @importFrom graphics plot grid points curve legend
#' @export
fit.plot <- function(fit, fit.type = 1, fit.legend = FALSE, l.size = 0.7) {
  if(!is.character(fit)) {
    stop("fit must be a character type object")
  }
  if(fit.type != 1 & fit.type != 2) {
    stop("fit.type must be equal to 1 or 2")
  }
  if(!is.numeric(l.size)) {
    stop("l.size must be numeric")
  }

  plot(NULL,  ylab = "Income shares", ylim = c(0, 1), xlim= c(0, 1),
    xlab = "Population shares", panel.first = grid(col = "gray78"))
  legendtext <- rep(NA, length(fit))
  for(i in 1:length(fit)){
    x <- eval(parse(text = fit[i]))
    if(fit.type == 1) {
      param <- t(as.matrix(x$nls.estimation[1, ]))
      colnames(param) <- colnames(x$nls.estimation)
      scale <- c("b", "mu")
      ind <- which(colnames(param) %in% scale)
      param[, colnames(param)[ind]] <- 10
    }
    if(fit.type == 2) {
      param <- t(as.matrix(x$gmm.estimation[1, ]))
    }
    legendtext[i] <- paste(c(x$distribution, " distribution"), collapse = "")
    shares <- c(0, x$grouped.data[1, ], 1)
    pop <- c(0, x$grouped.data[2, ], 1)
    points(pop, shares, pch = 20, type = "o", lty = 2)
    if(x$distribution == "GB2") {
      curve(lc.gb2(param, x), add = T, col = i+1)
    }
    if(x$distribution == "Beta 2") {
      curve(lc.gb2(c(1, param), x), add = T, col = i+1)
    }
    if(x$distribution == "Dagum") {
      curve(lc.gb2(c(param, 1), x), add = T, col = i+1)
    }
    if(x$distribution == "Singh-Maddala") {
      curve(lc.gb2(c(param[c(1, 2)], 1, param[3]), x), add = T, col = i+1)
    }
    if(x$distribution == "Fisk") {
      curve(lc.gb2(c(param[c(1, 2)], 1, 1), x), add = T, col = i+1)
    }
    if(x$distribution == "Lognormal") {
      curve(lc.ln(c(param[1], 1), x), add = T, col = i+1)
    }
  }
  if(fit.legend == TRUE) {
    legend("topleft", legend = legendtext, cex =  l.size, lty = 1,
      col = 2:(length(fit) + 1), ncol = 1)
  }
}

