#' Goodness of fit checks for binomial N-mixture models
#'
#' The package contains methods to compute overdispersion metrics, randomized quantile residuals,
#' and graphical diagnostics of model fit for binomial N-mixture models fitted using the \link[unmarked]{unmarked} package.
#' Details about the checks are given in Knape et al. (2018).
#'
#'
#' @docType package
#' @name nmixgof
#' @references Knape et al. 2018. Sensitivity of binomial N-mixture models to overdispersion: 
#' the importance of assessing model fit. Methods in Ecology and Evolution, 9:2102-2114. \doi{10.1111/2041-210X.13062}
#' @import unmarked
#' @useDynLib nmixgof
#' @importFrom Rcpp sourceCpp
#' @importFrom graphics abline
#' @importFrom stats dnbinom dpois pbinom plogis pnbinom ppois qnorm qqnorm runif
NULL


