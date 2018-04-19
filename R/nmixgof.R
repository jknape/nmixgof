#' Goodness of fit checks for binomial N-mixture models
#'
#' The package contains methods to compute overdispersion metrics, randomized quantile residuals,
#' and graphical diagnostics of model fit for binomial N-mixture models fitted using the \link[unmarked]{unmarked} package.
#'
#'
#' @docType package
#' @name nmixgof
#' @import unmarked
#' @useDynLib nmixgof
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline
#' @importFrom stats dnbinom dpois pbinom plogis pnbinom ppois qnorm qqnorm runif
NULL


