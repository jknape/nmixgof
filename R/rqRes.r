# Generic function for computing rq-residuals for count data from a CDF.
rqRes = function(y, pFun, ...) {
  b = pFun(y, ...)
  a = b*0
  a = pFun(y-1, ...)
  if (is.null(dim(y))) {
    res = numeric(length(y)) + NA
  } else {
    res = array(data= NA, dim = dim(y))
  }
  isnna = which(!is.na(y) & !is.na(a) & !is.na(b))
  res[isnna] = rqRes0(a[isnna], b[isnna])
  if (any(is.infinite(res))) {
    #browser()
    warning("Some residuals infinite.")
  }
  res
}

rqRes0 = function(a, b) {
  stopifnot(length(a) == length(b))
  qnorm(runif(length(a), a, b))
}

#' Plot residuals against fitted values
#'
#' Plots randomized-quantile residuals for binomial N-mixture models against fitted values.
#'
#' @param umFit An object from a model fitted using \link[unmarked]{pcount}.
#' @param type The type of randomized quantile residual to plot. One of 'margina', 'site-sum' or 'observation'.
#' @param residual The type of residual. Only 'rq' currently implemented.
#' @param ... Plot arguments.
#'
#' @export
#'
#' @examples
residfit = function(umFit, type = "marginal", residual = "rq", ...) {
  if (!inherits(umFit,"unmarkedFitPCount")) {
    stop("Argument needs to be a model fitted by the function pcount in the package unmarked.")
  }
  type = match.arg(type, c("marginal", "site-sum", "observation"))
  residual = match.arg(residual, "rq")
  if (identical(type, "marginal")) {
    rqr = rqResMar(umFit)
    fitted = fitted(umFit, K = umFit@K)
  } else if (identical(type, "site-sum")) {
    rqr = rqResSum(umFit)
    fitted = apply(fitted(umFit, K = umFit@K), 1, sum, na.rm = TRUE)
  } else {
    rqr = rqResObs(umFit)
    fitted = fitted(umFit, K = umFit@K)
  }
  plotArgs = list(fitted, rqr, ...)
  if (!("ylab" %in% names(plotArgs)))
      plotArgs = c(plotArgs, list(ylab = paste(type, "rq-residuals")))
  if (!("xlab" %in% names(plotArgs)))
    plotArgs = c(plotArgs, list(xlab = "fitted values"))
  do.call(plot, plotArgs)
}

#' Qq plot of randomized quantile residuals against standard normal quantiles
#'
#' @param umFit An object of class \link[unmarked]{unmarkedFitPCount} from a model fitted using \link[unmarked]{pcount}.
#' @param type The type of randomized quantile residual to plot. One of 'site-sum' or 'observation'.
#' @param main Plot label.
#' @param plotLine If true, the identity line is added to the plot.
#' @param ... Further arguments passed to qqnorm.
#'
#' @return A list with x and y coordinates of the qq plot, see \code[stats]{qqnorm}.
#' @export
#'
#' @examples
residqq =  function(umFit, type = "site-sum", main = "Residual qq plot", plotLine = TRUE, ...) {
    if (!inherits(umFit,"unmarkedFitPCount")) {
      stop("Argument needs to be a model fitted by the function pcount in the package unmarked.")
    }
    type = match.arg(type, c("marginal", "site-sum", "observation"))
    if (identical(type, "marginal")) {
      stop("Marginal qq-plot not implemented")
    } else if (identical(type, "site-sum")) {
      rqr = rqResSum(umFit)
    } else {
      rqr = rqResObs(umFit)
    }
    qqArgs = list(rqr ,main = main, ...)
    if (!("ylab" %in% names(qqArgs)))
      qqArgs = c(qqArgs , ylab = "residual quantile")
    qq = do.call(qqnorm, qqArgs)
    if (plotLine)
      abline(a=0, b=1)
    invisible(qq)
}

rqResMar = function(umFit) {
  fitval = fitted(umFit, K = umFit@K) # A bug in older versions of unmarked (fixed now) may cause incorrect fitted values for NB models unless K is supplied.
  y = umFit@data@y
  if (identical(umFit@mixture, "P")) {
    rqr = rqRes(y, pFun = ppois, lambda = fitval)
  } else if (identical(umFit@mixture, "NB")) {
    if (!identical(umFit@estimates["alpha"]@invlink, "exp"))
      stop("Unknown link function.")
    size = exp(umFit@estimates["alpha"]@estimates)
    rqr = rqRes(y, pFun = pnbinom, mu = fitval, size = size)
  } else if (identical(umFit@mixture, "ZIP")) {
    if (!identical(umFit@estimates["psi"]@invlink, "logistic"))
      stop("Unknown link function.")
    pZIP = function(y, lambda, psi) {
      psi*(y>=0) + (1-psi) * ppois(y, lambda/(1-psi))
    }
    psi = plogis(umFit@estimates["psi"]@estimates)
    rqr = rqRes(y, pFun = pZIP, lambda = fitval, psi = psi)
  } else {stop("Mixture not recognized.")}
  rqr
}

rqResObs = function(umFit) {
  rN = integer(nrow(umFit@data@y)) + NA
  if (length(umFit@sitesRemoved) > 0)
    rN[-umFit@sitesRemoved] = apply(unmarked::ranef(umFit)@post[,,1], 1 , sample, x =umFit@K + 1, size = 1, replace = FALSE) - 1
  else
    rN = apply(unmarked::ranef(umFit)@post[,,1], 1 , sample, x =umFit@K + 1, size = 1, replace = FALSE) - 1
  p = getP(umFit)
  res = rqRes(umFit@data@y, pFun = pbinom, size = kronecker(rN, t(rep(1, ncol(p)))), prob=p)
  if (any(is.infinite(res))) {
    #browser()
    warning(paste(sum(is.infinite(res)), " residuals infinite."))
  }
  res
}


rqResSum = function(umFit) {
  lam = predict(umFit, type="state")[,1]
  p = getP(umFit)
  res = numeric(nrow(umFit@data@y)) + NA
  if (length(umFit@sitesRemoved) > 0)
    y = umFit@data@y[-umFit@sitesRemoved,]
  else
    y = umFit@data@y
  cumProb = matrix(0, nrow = nrow(y), ncol = 2)
  dfun = switch(umFit@mixture,
                P = function(N) {dpois(N, lam)},
                NB = function(N) {dnbinom(N, mu=lam, size=exp(coef(umFit, type="alpha")))},
                ZIP = function(N) {psi = plogis(coef(umFit, type="psi"))
                                   (1-psi)*dpois(N, lam/(1-psi)) + psi*(N==0)}
         )
  for (N in 0:umFit@K) {
    cumProb = cumProb + kronecker(dfun(N), t(c(1,1))) * pbinsum(y, rep(N, nrow(y)), p)[,2:3]
  }
  if (length(umFit@sitesRemoved) > 0)
    res[-umFit@sitesRemoved] = rqRes0(cumProb[,1], cumProb[,2])
  else
    res = rqRes0(cumProb[,1], cumProb[,2])
  if (any(is.infinite(res))) {
    #browser()
    warning(paste(sum(is.infinite(res)), " residuals infinite."))
  }
  res
}

plotgof = function(umFit, QQ = c("site-sum", "obs", "marginal"), covariate ) {
  QQ = match.arg(QQ, c("site-sum", "obs", "marginal", ""), several.ok = TRUE)
  rqrS = rqResSum(umFit)
  if ("site-sum" %in% QQ) {
    qqnorm(rqrS, main = "QQ-plot for site-sum rq-residuals")
    abline(a=0, b=1)
  }
  rqrO = rqResObs(umFit)
  if ("obs" %in% QQ) {
    qqnorm(rqrO, main = "QQ-plot for observation rq-residuals")
    abline(a=0, b=1)
  }
  data = unmarked::getData(umFit)
  if (length(umFit@sitesRemoved) > 0)
    siteCovs = data@siteCovs[-umFit@sitesRemoved, ]
  else
    siteCovs = data@siteCovs
  for (i in 1:ncol(siteCovs)) {
    resCovPlot(siteCovs[, i], rqrS, xlab = colnames(siteCovs)[i], ylab = "site sum rq-residual")
  }
  obsCovs = data@obsCovs
  rqrM = rqResMarginal(umFit)
  for (i in 1:ncol(obsCovs)) {
    resCovPlot(obsCovs[, i], as.vector(rqrM), xlab = colnames(obsCovs)[i], ylab = "marginal rq-residual")
  }

  for (i in 1:ncol(obsCovs)) {
    resCovPlot(obsCovs[, i], as.vector(rqrO), xlab = colnames(obsCovs)[i], ylab = "obs rq-residual")
  }
}


#' Overdispersion metrics for binomial N-mixture models.
#'
#' Computes various types of overdispersion metrics for binomial N-mixture models.
#' @param umFit An object of class \link[unmarked]{unmarkedFitPCount} from a model fitted using \link[unmarked]{pcount}.
#' @param type The type of metric to compute, one of 'marginal', 'site-sum' or 'observation'.
#'
#' @return An estimate of overdispersion relative to the fitted model.
#' @export
#'
#' @examples
#' library(unmarked)
#' data(mallard)
#' fm.mallard <- pcount(~ 1 ~ 1, unmarkedFramePCount(y = mallard.y), K=100)
#' chat(fm.mallard, "m")
#' chat(fm.mallard, "s")
#' chat(fm.mallard, "o")
chat = function(umFit, type = "marginal") {
  if (!inherits(umFit,"unmarkedFitPCount")) {
    stop("Argument needs to be a model fitted by the function pcount in the package unmarked.")
  }
  type = match.arg(type, c("marginal", "site-sum", "observation"))
  if (identical(type, "marginal"))
    chat = chatM(umFit)
  else if (identical(type, "site-sum"))
    chat = chatS(umFit)
  else
    chat = chatO(umFit)
  chat
}

chatS = function(umFit) {
  p = getP(umFit)
  lam = predict(umFit, "state")[,1]
  odPar = switch(umFit@mixture, P = NA, ZIP = coef(umFit)["psi(psi)"], NB  = coef(umFit)["alpha(alpha)"])
  nVar = switch(umFit@mixture, P = lam, ZIP = lam*(1 + lam * plogis(odPar)/(1-plogis(odPar))), NB = lam + lam^2 * exp(-odPar))
  if (length(umFit@sitesRemoved) > 0) {
    y = umFit@data@y[-umFit@sitesRemoved,]
    expected = fitted(umFit, K = umFit@K)[-umFit@sitesRemoved,]
  } else {
    y = umFit@data@y
    expected = fitted(umFit, K = umFit@K)
  }
  naMat = (is.finite(y + expected))
  naMat[which(naMat != 1)] = NA
  obs.site = apply(y * naMat, 1, sum, na.rm = TRUE)
  exp.site = apply(expected * naMat, 1, sum, na.rm = TRUE)
  var.site = lam * apply(p * (1 - p) * naMat, 1, sum, na.rm = TRUE) +
    nVar * apply(p * naMat, 1, function(row) {sum(outer(row, row), na.rm = TRUE)})
  chi2 = sum((obs.site - exp.site)^2/var.site)
  chi2/(nrow(y) - length(coef(umFit)))
}

chatO = function(umFit) {
  p = getP(umFit)
  if (length(umFit@sitesRemoved) > 0) {
    y = umFit@data@y[-umFit@sitesRemoved,]
  } else {
    y = umFit@data@y
  }
  rN = apply(unmarked::ranef(umFit)@post[,,1], 1 , sample, x =umFit@K + 1, size = 1, replace = FALSE) - 1
  naMat = (is.finite(y + p))
  naMat[which(naMat != 1)] = NA
  obs.site = apply(y * naMat, 1, sum, na.rm = TRUE)
  exp.site = apply(p * naMat, 1, sum, na.rm = TRUE) *rN
  var.site = apply(p * (1 - p) * naMat, 1, sum, na.rm = TRUE) * rN
  posY = which(obs.site>0)
  chi2 = sum((obs.site[posY] - exp.site[posY])^2/var.site[posY])
  chi2/(length(posY) - length(coef(umFit)))
}

chatM = function(umFit) {
  odPar = switch(umFit@mixture, P = NA, ZIP = coef(umFit)["psi(psi)"], NB  = coef(umFit)["alpha(alpha)"])
  if (length(umFit@sitesRemoved) > 0) {
    y = umFit@data@y[-umFit@sitesRemoved,]
    expected = fitted(umFit, K = umFit@K)[-umFit@sitesRemoved,]
  } else {
    y = umFit@data@y
    expected = fitted(umFit, K = umFit@K)
  }
  ## For ZIP, zero inflation has been included in expected.
  variance = switch(umFit@mixture, P = expected, ZIP = expected*(1 + expected * plogis(odPar)/(1-plogis(odPar))), NB = expected + expected^2 * exp(-odPar))
  res = (y-expected)/sqrt(variance)
  chi2 = sum(res^2, na.rm = TRUE)
  chi2/(sum(!is.na(res)) - length(coef(umFit)))
}


#' Plot residuals against covariates
#'
#'
#' Site-sum randomized quantile residuals are used for site covariates while marginal residuals are used for observation covariates.
#' The same random residual draws are reused for different covariates.
#'
#' @param umFit An object of class \link[unmarked]{unmarkedFitPCount} from a model fitted using \link[unmarked]{pcount}.
#' @param ...
#'
#' @export
#'
#' @examples
residcov = function(umFit, ...) {
  data = getData(umFit)
  if (!is.null(data@siteCovs)) {
    rqS = rqResSum(umFit)
    for (i in 1:ncol(data@siteCovs)) {
      plot(data@siteCovs[,i], rqS, xlab = colnames(data@siteCovs)[i],ylab = "site-sum rq resiudal", ...)
    }
  }
  if (!is.null(data@obsCovs)) {
    rqM = rqResMar(umFit)
    for (i in 1:ncol(data@obsCovs)) {
      obsCov = matrix(data@obsCovs[,i], nrow = nrow(data@y), ncol = ncol(data@y),byrow = TRUE)
      plot(data@obsCovs[,i], rqM, xlab = colnames(data@obsCovs)[i],ylab = "marginal rq resiudal", ...)
    }
  }
}

