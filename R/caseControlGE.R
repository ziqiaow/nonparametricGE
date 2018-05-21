#' Semiparametric Gene-Environment Interactions in Case-Control Studies
#'
#' This package was written to perform simulations for the follow up to Stalder et. al. (2017).
#'
#' \code{\link{simfun}} Simulates then maximizes the retrospective pseudolikelihood of case-control data under the
#' assumption of G-E independence.  It returns estimates, asymptotic standard error, and
#' (optionally) bootstrapped standard error for several versions of the estimator.
#'
#' \code{\link{simulate_complex}} simulates data to be analyzed by \code{simfun},
#' logistic regression, or other models.  This function was used to simulate data for
#' Stalder et. al. (2017).  The user can specify up to three types of genetic variables
#' (SNPs with additive effects under HWE, normally distributed polygenic risk scores,
#' and gamma distributed polygenic risk scores), each of which can be multivariate.
#' Two types of environmental variables (binary and normal) can also be potentially multivariate.
#'
#' \code{\link{combo_asymp}} is called by \code{\link{simfun}} to calculate estimates and asymptotic
#' SE.
#'
#' \code{\link{combo_boot}} is called by \code{\link{simfun}} to calculate bootstrap estimates.  It
#' is called after \code{\link{combo_asymp}}
#'
#' @useDynLib caseControlGE
#' @importFrom Rcpp evalCpp
#' @importFrom utils head
#' @importFrom stats binomial coef cor cov glm lm model.matrix plogis qlogis qnorm rbinom rgamma rnorm rt runif sd
#' @importFrom stats model.frame setNames
"_PACKAGE"
