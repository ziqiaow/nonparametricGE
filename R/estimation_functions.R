##### SPMLE functions, Alex Asher, 28 May, 2018

## Wrapper to calculate negative loglikelihood and gradient.  Report likelihood, store gradient
lik_fn = function(Omega, D, G, E, pi1) {
  lik_grad = neglikgrad(Omega=Omega, D=D, G=G, E=E, pi1=pi1) # neglikgrad is c++ code that calculates the negative loglikelihood and gradient
  assign(x="grad_val", value=lik_grad[["gradient"]], pos=parent.frame())  # Store the gradient in the calling environment to retrieve it when ucminf tries to calculate the gradient
  return(lik_grad[["objective"]])
}


## Report the gradient
grad_fn = function(Omega, D, G, E, pi1) return(get(x="grad_val", pos=parent.frame())) # return the previously stored gradient


## When ucminf works, it works brilliantly (typically more than twice as fast
## as the next-fastest algorithm).  But it has a nasty habit of declaring
## convergence before actually converging.
## max_grad_tol = maximum allowable gradient
## num_retries = number of times to retry optimization
##
## Omega_start is the starting value for optimization and invhessian.lt is the lower triangle of the inverse hessian

## Calculate the SPMLE (optionally precondition with hessian)
maximize_spmle = function(Omega_start, D, G, E, pi1, control=list()) {
  ## Convert G, and E to matrices (if not already)
  G = as.matrix(G)
  E = as.matrix(E)

  ## Set control parameters for SPMLE estimation
  con = list(trace=0, max_grad_tol=0.001, num_retries=2)
  con[(names(control))] = control

  ## Set control parameters for ucminf using values from con where appropriate
  ucminf_con = con[names(con) %in% c("trace", "grtol", "xtol", "stepmax", "maxeval", "grad", "gradstep", "invhessian.lt")]

  ## Optimize with ucminf
  spmle_max = ucminf::ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)

  ## Save the retry number and starting values
  spmle_max$start_vals = Omega_start
  spmle_max$retry = 0

  ## Check for convergance and rerun optimization if necessary
  if(spmle_max$convergence<1 || !is.finite(spmle_max$info[1]) || spmle_max$info[1]>con$max_grad_tol) {
    ## Scale starting values by 1/SD
    scale_vec = 1/c(1, apply(model.matrix(~0+G*E), 2, sd))

    ## Loop through different starting values
    for(j in seq_len(con$num_retries)) {
      if(con$trace>-1) {cat("ucminf retry", j, "of", con$num_retries, "\n")}

      ## If failure happened when preconditioning with the hessian, try without (same start values the first time)
      if(!is.null(ucminf_con$invhessian.lt)) {
        ucminf_con$invhessian.lt = NULL
        new_start = Omega_start
      } else {  # otherwise use scaled normal startvals (preserve names from Omega_start)
        new_start = setNames(rnorm(n=length(Omega_start), sd=scale_vec), names(Omega_start))
      }

      ## Rerun optimization with new values
      spmle_max = ucminf::ucminf(par=new_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)

      ## Save the retry number and starting values
      spmle_max$start_vals = new_start
      spmle_max$retry = j

      ## Break out of the loop once we have convergence
      if(is.finite(spmle_max$info[1]) & spmle_max$info[1]<con$max_grad_tol) {break}
    }

    ## Test again: if still no convergance, throw an error
    if(!is.finite(spmle_max$info[1]) | spmle_max$info[1]>con$max_grad_tol) {stop("ucminf failed to converge. Try increasing control$num_retries or rescaling G and E.")}
  }
  return(spmle_max)
}


#' Semiparametric Maximum Pseudolieklihood Estimator for Case-Control Studies Under G-E Independence.
#'
#' \code{spmle} maximizes the retrospective pseudolikelihood of case-control data under the assumption
#' of G-E independence in the underlying population.
#' The marginal distributions of G and E are treated nonparametrically.
#'
#' This function applies the method of Stalder et. al. (2017) to maximize the
#' retrospective pseudolikelihood of case-control data under the assumption of G-E independence.
#' It currently supports the model with G and E main effects, and a multiplicative G*E interaction.
#'
#' The \code{control} argument is a list that controls the behavior of the optimization algorithm
#' \code{\link[ucminf]{ucminf}} from the \pkg{ucminf} package.  When \code{ucminf} works,
#' it works brilliantly (typically more than twice as fast as the next-fastest algorithm).
#' But it has a nasty habit of declaring convergence before actually converging.
#' To address this, \code{spmle} checks the maximum gradient at "convergence", and can rerun the optimization
#' using different starting values.  The \code{control} argument can supply any of the following components:
#' \describe{
#'   \item{\code{max_grad_tol}}{maximum allowable gradient at convergence.  \code{spmle} does not
#'     consider the optimization to have converged if the maximum gradient \code{> max_grad_tol}.
#'     Default \code{max_grad_tol = 0.001}.}
#'   \item{\code{num_retries}}{number of times to retry optimization.  An error is produced if
#'     the optimization has not converged after \code{num_retries}.  Different starting values
#'     are used for each retry.  Default \code{num_retries = 2}.}
#'   \item{\code{use_hess}}{a logical value instructing \code{spmle} to use the analytic hessian
#'     to precondition the optimization.  This brings significant speed benefits, and is one reason
#'     \code{ucminf} is so fast.  For unknown reasons, preconditioning causes computers with
#'     certain Intel CPUs to prematurely terminate iterating.  By default, \code{use_hess = TRUE},
#'     but if you notice that \code{ucminf} never converges during the first attempt, try setting \code{use_hess = FALSE}.}
#'   \item{\code{trace}}{a scalar or logical value that is used by both \code{spmle} and \code{ucminf}
#'     to control the printing of detailed tracing information.
#'     If TRUE or >0, details of each \code{ucminf} iteration are printed.
#'     If FALSE or 0, \code{ucminf} iteration details are suppressed but \code{spmle} still
#'     prints optimization retries.  If \code{trace < 0} nothing is printed.  Default \code{trace = 0}.}
#'   \item{additional control parameters}{not used by \code{spmle}, but are passed to \code{\link[ucminf]{ucminf}}.
#'     Note that the \code{ucminf} algorithm has four stopping criteria, and \code{ucminf} will
#'     declare convergence if any one of them has been met.  The \code{ucminf} control parameter
#'     "\code{grtol}" controls \code{ucminf}'s gradient stopping criterion, which defaults to
#'     \code{1e-6}.  \code{grtol} should not be set larger than the \code{spmle} control parameter \code{max_grad_tol}.}
#' }
#'
#' @section References:
#' Stalder, O., Asher, A., Liang, L., Carroll, R. J., Ma, Y., and Chatterjee, N. (2017).
#' \emph{Semi-parametric analysis of complex polygenic gene-environment interactions in case-control studies.}
#' Biometrika, 104, 801–812.
#'
#' @param D a binary vector of disease status (1=case, 0=control).
#' @param G a vector or matrix (if multivariate) containing genetic data. Can be continuous, discrete, or a combination.
#' @param E a vector or matrix (if multivariate) containing environmental data. Can be continuous, discrete, or a combination.
#' @param pi1 the population disease rate, a scalar in [0, 1) or the string "rare".
#'   Using \code{pi1=0} is the rare disease approximation.
#' @param data an optional list, environment, or object coercible to a data frame by
#'   \code{\link[base]{as.data.frame}} containing the variables in the model.  If not found
#'   in data, the variables are taken from the environment from which \code{spmleCombo} is called.
#' @param control a list of control parameters that allow the user to control the optimization algorithm.  See 'Details'.
#' @param swap a logical scalar rarely of interest to the end user.  Dependence on the distributions of G and E
#'   are removed using different methods; this switch swaps them to produce a symmetric estimator with identical
#'   properties to the SPMLE.  Default \code{FALSE}.
#' @param startvals an optional numeric vector of coefficient starting values for optimization.  Usually left blank,
#'   in which case logistic regression estimates are used as starting values.
#' @return an object of class \code{"spmle"}.  The function \code{summary} (i.e., \code{summary.spmle})
#'   can be used to obtain or print a summary of the results and the function \code{anova}
#'   (i.e., anova.spmle) to produce an analysis of variance table using likelihood-ratio tests.
#'   Note that \code{anova.spmle} may be used to compare one \code{"spmle"} object
#'   to another, because the loglikelihood reported by \code{logLik(spmle-obj)} is
#'   accurate up to an additive constant; however it should not be used to compare
#'   an \code{spmle-object} to a model fit by a different method.
#'
#'   \code{\link{predict.spmle}}, the \code{predict} method for S3 class \code{"spmle"}
#'   can predict the expected response (on logistic or probability scales), compute
#'   confidence intervals for the expected response, and provide standard errors.
#'
#'   The generic accessor functions \code{coefficients}, \code{fitted.values}
#'   and \code{residuals} can be used to extract various useful features of the value returned by spmle.
#'
#'   An object of class \code{"spmle"} is a list containing at least the following components:
#' \describe{
#'   \item{\code{coefficients}}{a named vector of coefficients}
#'   \item{\code{pi1}}{the value of pi1 used during the analysis}
#'   \item{\code{SE}}{standard error estimate of coefficients}
#'   \item{\code{cov}}{covariance matrix of coefficients}
#'   \item{\code{glm_fit}}{a logistic regression model fit using the same model as \code{spmle}}
#'   \item{\code{call}}{the matched call}
#'   \item{\code{formula}}{the formula supplied}
#'   \item{\code{data}}{the \code{data argument}}
#'   \item{\code{model}}{the model frame}
#'   \item{\code{terms}}{the \code{terms} object used}
#'   \item{\code{linear.predictors}}{the linear fit on the logistic link scale}
#'   \item{\code{fitted.values}}{the fitted values on the probability scale}
#'   \item{\code{residuals}}{the Pearson residuals}
#'   \item{\code{null.deviance}}{the deviance for the null model.  Deviance = \code{-2*logLik}.}
#'   \item{\code{df.residual}}{the residual degrees of freedom}
#'   \item{\code{df.null}}{the residual degrees of freedom for the null model}
#'   \item{\code{rank}}{the numeric rank of the fitted linear model (i.e. df_model:
#'     the number of parameters estimated)}
#'   \item{\code{nobs}}{number of observations}
#'   \item{\code{ncase}}{number of cases}
#'   \item{\code{ncontrol}}{number of controls}
#' }
#'
#' \code{spmle} objects created by \code{spmle()} additionally have components \code{logLik}
#' (log pseudolikelihood), \code{deviance} (-2 * log pseudolikelihood), \code{aic}, \code{bic},
#' \code{ucminf} (optimization output), and matrices \code{H_inv}, \code{Sigma}, \code{zeta0},
#' and \code{zeta1}, which are used in calculating the asymptotic estimate of standard error.
#'
#' @seealso \code{\link{spmleCombo}} for a slower but more precise estimator, \code{\link{simulateCC}} to simulate data
#' @examples
#' # Simulation from Table 1 in Stalder et. al. (2017)
#' set.seed(2018)
#' dat = simulateCC(ncase=500, ncontrol=500, beta0=-4.165,
#'                  betaG_SNP=c(log(1.2), log(1.2), 0, log(1.2), 0),
#'                  betaE_bin=log(1.5),
#'                  betaGE_SNP_bin=c(log(1.3), 0, 0, log(1.3), 0),
#'                  MAF=c(0.1, 0.3, 0.3, 0.3, 0.1),
#'                  SNP_cor=0.7, E_bin_freq=0.5)
#'
#' # SPMLE with known population disease rate of 0.03
#' spmle(D=D, G=G, E=E, pi1=0.03, data=dat)
#'
#' # Simulation with a single SNP and a single binary environmental variable.
#' # True population disease rate in this simulation is 0.03.
#' # This simulation scenario was used in the Supplementary Material of Stalder et. al. (2017)
#' # to compare performance against the less flexible method of Chatterjee and Carroll (2005),
#' # which is available as the function as snp.logistic in the Bioconductor package CGEN.
#' dat2 = simulateCC(ncase=100, ncontrol=100, beta0=-3.77,
#'                   betaG_SNP=log(1.2), betaE_bin=log(1.5),
#'                   betaGE_SNP_bin=log(1.3), MAF=0.1,
#'                   E_bin_freq=0.5)
#'
#' # SPMLE using the rare disease assumption, optimization tracing,
#' # and no hessian preconditioning.
#' spmle(D=D, G=G, E=E, pi1=0, data=dat2, control=list(trace=0, use_hess=FALSE))
#' @export
spmle = function(D, G, E, pi1, data, control=list(), swap=FALSE, startvals){
  ## Store the function call
  cl = match.call()

  ## Get argument names for D, G, and E
  Dname = substitute(D)
  Gname = substitute(G)
  Ename = substitute(E)

  ## If pi1 (partially) matches "rare", set pi1=0
  if(pmatch(x=substr(x=tolower(pi1), start=1, stop=4), table="rare", nomatch=FALSE)) {pi1=0}
  stopifnot(pi1>=0, pi1<1)  # throw an error if pi1 not in [0,1)

  ## Store the formula with user-provided variable names.  For consistency with estimators that accept formulas (like lm),
  ## set the formula environment to the calling environment (as if formula had been an argument).
  formula = formula(paste(as.character(as.expression(Dname)), "~", as.character(as.expression(Gname)), "*", as.character(as.expression(Ename))))
  attr(formula, ".Environment") = parent.frame()

  ## If no data.frame was supplied, set data to the calling environment
  if(missing(data)) {
    data = environment(formula)
  } else {                          # evaluate D, G, and E in data, if appropriate
    if(class(data)=="matrix") {data = as.data.frame(data)}
    D = with(data=data, eval(Dname))
    G = with(data=data, eval(Gname))
    E = with(data=data, eval(Ename))
  }

  ## Expand factor variables into dummies, if present
  if(class(G)=="factor") {G=model.matrix(~G)[,-1]}
  if(class(E)=="factor") {E=model.matrix(~E)[,-1]}

  ## Create model.frame (response, no interactions) and model.matrix (all predictors, including interactions)
  model_frame = model.frame(formula=formula, data=data)
  model_matrix = model.matrix(formula, model_frame)

  ## Fit the model with logistic regression
  logistic_fit = glm(formula, family=binomial(link='logit'), data=data)

  ## If starting values weren't provided, use logistic regression estimates
  if(missing(startvals)) {
    Omega_start = coef(logistic_fit)
  } else {
    Omega_start = startvals
  }

  ## If user-provided startval lacked names, add them
  if(is.null(names(Omega_start))) {names(Omega_start) = colnames(model_matrix)}

  ## Set control parameters
  con = list(trace=0, use_hess=TRUE, max_grad_tol=0.001, num_retries=2)
  con[(names(control))] = control

  ## Remove missing values and convert G and E to matrices
  complete = complete.cases(D, G, E)
  D = D[complete]
  G = as.matrix(as.matrix(G)[complete,])
  E = as.matrix(as.matrix(E)[complete,])

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  G = as.matrix(G)
  E = as.matrix(E)

  ## Model degrees of freedom
  df_model = length(Omega_start)

  ## If we're swapping G & E, make the change now
  if(swap==TRUE) {
    nG = NCOL(G)
    nE = NCOL(E)
    swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))
    reverse_swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
    temp = G
    G = E
    E = temp
  } else {
    swap_order = reverse_swap_order = seq_along(Omega_start)
  }

  ## Precondition with hessian if requested
  if(con$use_hess==TRUE) {con$invhessian.lt = solve(neghess(Omega=Omega_start, D=D, G=G, E=E, pi1=pi1))[lower.tri(diag(df_model),diag=TRUE)]}

  ## Calculate SPMLE for a given disease rate
  spmle_max = maximize_spmle(Omega_start=Omega_start[swap_order], D=D, G=G, E=E, pi1=pi1, control=con)

  ## Asymptotic SE for known disease rate
  hess_zeta = hesszeta(Omega=spmle_max$par, D=D, G=G, E=E, pi1=pi1)
  Sigma = ((ncontrol-1)*cov(hess_zeta$zeta0) + (ncase-1)*cov(hess_zeta$zeta1))/n  # (ncontrol-1) to correct denominator because cov uses the "sample covariance matrix" denominator of n-1
  H_inv = -n*solve(hess_zeta$hessian)    # = (-hessian/n)^-1 = (Gamma_1 - Gamma_2)^-1
  Lambda = H_inv %*% Sigma %*% t(H_inv) / n  # covar matrix of sqrt(n) * OmegaHat
  SE_asy = sqrt(diag(Lambda))

  ## Calculate predictions
  linear_predictors = model_matrix %*% spmle_max$par  # logistic scale
  fitted_values = plogis(q=linear_predictors)         # probability scale

  ## Loglikelihood, AIC, and BIC
  loglik = -spmle_max$value
  attr(loglik, "nobs") = n
  attr(loglik, "df") = df_model
  class(loglik) = "logLik"
  AIC = 2*df_model - 2*loglik
  BIC = log(n) * df_model - 2*loglik

  ## Pearson residuals, total deviance, and df resid & null
  pearson_resid = (D-fitted_values) / sqrt(fitted_values*(1-fitted_values))
  total_deviance = -2*loglik
  df_resid = n - df_model
  null_deviance = sum(-2*log(abs(1-D-mean(D))))
  df_null = n - 1

  ## Compile results into a list.  Use glm-object naming conventions.  If we swapped G & E, change back now
  spmle_est = list(coefficients      = spmle_max$par[reverse_swap_order],
                   pi1               = pi1,
                   SE                = SE_asy[reverse_swap_order],
                   cov               = Lambda[reverse_swap_order,reverse_swap_order],
                   H_inv             = H_inv[reverse_swap_order,reverse_swap_order],
                   Sigma             = Sigma[reverse_swap_order,reverse_swap_order],
                   zeta0             = hess_zeta$zeta0[, reverse_swap_order],
                   zeta1             = hess_zeta$zeta1[, reverse_swap_order],
                   ucminf            = spmle_max,
                   glm_fit           = logistic_fit,
                   call              = cl,
                   formula           = formula,
                   data              = data,
                   model             = model_frame,
                   terms             = attr(model_frame, "terms"),
                   linear.predictors = linear_predictors,
                   fitted.values     = fitted_values,
                   residuals         = pearson_resid,
                   deviance          = total_deviance,
                   null.deviance     = null_deviance,
                   df.residual       = df_resid,
                   df.null           = df_null,
                   rank              = df_model,
                   aic               = AIC,
                   bic               = BIC,
                   logLik            = loglik,
                   nobs              = n,
                   ncase             = ncase,
                   ncontrol          = ncontrol)

  ## Return an object of class spmle
  class(spmle_est) = c("spmle")
  return(spmle_est)
}


#' Improved Semiparametric Estimator for Case-Control Studies Under G-E Independence.
#'
#' \code{spmleCombo} estimates Gene-Environment interactions in case-control data under the assumption
#' of G-E independence in the underlying population.  This is an improved version of
#' \code{\link{spmle}} that improves estimation accuracy at the cost of
#' computing time.
#'
#' This function calculates the Symmetric Combination Estimator of Wang et. al. (2018) to combine two symmetric
#' retrospective pseudolikelihood estimators of case-control data under the assumption of G-E independence.
#' The marginal distributions of G and E are treated nonparametrically.  This is done by
#' calling \code{spmle} twice (once with \code{swap=TRUE}) to generate two symmetric
#' estimates, which are then combined using a GLS approach.
#' It currently supports the model with G and E main effects, and a multiplicative G*E interaction.
#'
#' The \code{control} argument is a list that controls the behavior of the optimization algorithm
#' \code{\link[ucminf]{ucminf}} from the \pkg{ucminf} package.  When \code{ucminf} works,
#' it works brilliantly (typically more than twice as fast as the next-fastest algorithm).
#' But it has a nasty habit of declaring convergence before actually converging.
#' To address this, \code{spmleCombo} checks the maximum gradient at "convergence", and can rerun the optimization
#' using different starting values.  The \code{control} argument can supply any of the following components:
#' \describe{
#'   \item{\code{max_grad_tol}}{maximum allowable gradient at convergence.  \code{spmleCombo} does not
#'     consider the optimization to have converged if the maximum gradient \code{> max_grad_tol}.
#'     Default \code{max_grad_tol = 0.001}.}
#'   \item{\code{num_retries}}{number of times to retry optimization.  An error is produced if
#'     the optimization has not converged after \code{num_retries}.  Different starting values
#'     are used for each retry.  Default \code{num_retries = 2}.}
#'   \item{\code{use_hess}}{a logical value instructing \code{spmleCombo} to use the analytic hessian
#'     to precondition the optimization.  This brings significant speed benefits, and is one reason
#'     \code{ucminf} is so fast.  For unknown reasons, preconditioning causes computers with
#'     certain Intel CPUs to prematurely terminate iterating.  By default, \code{use_hess = TRUE},
#'     but if you notice that \code{ucminf} never converges during the first attempt, try setting \code{use_hess = FALSE}.}
#'   \item{\code{trace}}{a scalar or logical value that is used by both \code{spmleCombo} and \code{ucminf}
#'     to control the printing of detailed tracing information.
#'     If TRUE or >0, details of each \code{ucminf} iteration are printed.
#'     If FALSE or 0, \code{ucminf} iteration details are suppressed but \code{spmleCombo} still
#'     prints optimization retries.  If \code{trace < 0} nothing is printed.  Default \code{trace = 0}.}
#'   \item{additional control parameters}{not used by \code{spmleCombo}, but are passed to \code{\link[ucminf]{ucminf}}.
#'     Note that the \code{ucminf} algorithm has four stopping criteria, and \code{ucminf} will
#'     declare convergence if any one of them has been met.  The \code{ucminf} control parameter
#'     "\code{grtol}" controls \code{ucminf}'s gradient stopping criterion, which defaults to
#'     \code{1e-6}.  \code{grtol} should not be set larger than the \code{spmleCombo} control parameter \code{max_grad_tol}.}
#' }
#'
#' @section References:
#' Stalder, O., Asher, A., Liang, L., Carroll, R. J., Ma, Y., and Chatterjee, N. (2017).
#' \emph{Semi-parametric analysis of complex polygenic gene-environment interactions in case-control studies.}
#' Biometrika, 104, 801–812.
#'
#' Wang, T., Asher, A., Carroll, R. J. (2018).
#' \emph{Improved Semiparametric Analysis of Polygenic Gene-Environment Interactions in Case-Control Studies}
#' Unpublished.
#'
#' @param D a binary vector of disease status (1=case, 0=control).
#' @param G a vector or matrix (if multivariate) containing genetic data. Can be continuous, discrete, or a combination.
#' @param E a vector or matrix (if multivariate) containing environmental data. Can be continuous, discrete, or a combination.
#' @param pi1 the population disease rate, a scalar in [0, 1) or the string "rare".
#'   Using \code{pi1=0} is the rare disease approximation.
#' @param data an optional list, environment, or object coercible to a data frame by
#'   \code{\link[base]{as.data.frame}} containing the variables in the model.  If not found
#'   in data, the variables are taken from the environment from which \code{spmleCombo} is called.
#' @param nboot an integer: the number of bootstraps to use when estimating the SE of the
#'   Symmetric Combination Estimator.  Setting \code{nboot=0} disables the bootstrap
#'   and uses the asymptotic standard error estimate (not recommended because asymptotic
#'   estimates of the SE can provide poor coverage - setting \code{nboot=0} will
#'   trigger a warning).  Default \code{50}.
#' @param ncores an integer: the number of cpu cores to use when parallelizing bootstraps.
#'   Default \code{1} executes bootstrap sequentially.
#' @param control a list of control parameters that allow the user to control the optimization algorithm.  See 'Details'.
#' @param startvals an optional numeric vector of coefficient starting values for optimization.  Usually left blank,
#'   in which case logistic regression estimates are used as starting values.
#' @return an object of class \code{"spmle"}.  The function \code{summary} (i.e., \code{summary.spmle})
#'   can be used to obtain or print a summary of the results.  The Symmetric Combination Estimator
#'   is not a maximum (pseudo)likelihood estimator like \code{\link{spmle}}; it is the
#'   optimal combination of two such estimators.  As such, it has no associated loglikelihood
#'   and the function \code{anova.spmle} cannot be used to compare models fit using
#'   \code{spmleCombo}.
#'
#'   \code{\link{predict.spmle}}, the \code{predict} method for S3 class \code{"spmle"}
#'   can predict the expected response (on logistic or probability scales), compute
#'   confidence intervals for the expected response, and provide standard errors.
#'
#'   The generic accessor functions \code{coefficients}, \code{fitted.values}
#'   and \code{residuals} can be used to extract various useful features of the value returned by spmle.
#'
#'   An object of class \code{"spmle"} is a list containing at least the following components:
#' \describe{
#'   \item{\code{coefficients}}{a named vector of coefficients}
#'   \item{\code{pi1}}{the value of pi1 used during the analysis}
#'   \item{\code{SE}}{standard error estimate of coefficients}
#'   \item{\code{cov}}{covariance matrix of coefficients}
#'   \item{\code{glm_fit}}{a logistic regression model fit using the same model as \code{spmleCombo}}
#'   \item{\code{call}}{the matched call}
#'   \item{\code{formula}}{the formula supplied}
#'   \item{\code{data}}{the \code{data argument}}
#'   \item{\code{model}}{the model frame}
#'   \item{\code{terms}}{the \code{terms} object used}
#'   \item{\code{linear.predictors}}{the linear fit on the logistic link scale}
#'   \item{\code{fitted.values}}{the fitted values on the probability scale}
#'   \item{\code{residuals}}{the Pearson residuals}
#'   \item{\code{null.deviance}}{the deviance for the null model}
#'   \item{\code{df.residual}}{the residual degrees of freedom}
#'   \item{\code{df.null}}{the residual degrees of freedom for the null model}
#'   \item{\code{rank}}{the numeric rank of the fitted linear model (i.e. df_model:
#'     the number of parameters estimated)}
#'   \item{\code{nobs}}{number of observations}
#'   \item{\code{ncase}}{number of cases}
#'   \item{\code{ncontrol}}{number of controls}
#' }
#'
#' \code{spmle} objects created by \code{spmleCombo()} additionally have components \code{spmle_E}
#' (model from \code{spmle} that profiled out the distribution of E), \code{spmle_G}
#' (model from \code{spmle} that profiled out the distribution of G with \code{swap=TRUE}),
#' and \code{bootstraps} (matrix of bootstrapped parameter estimates, if \code{nboot > 0}).
#'
#' @seealso \code{\link{spmle}}, \code{\link{simulateCC}} to simulate data
#' @examples
#' # Simulation from Table 1 in Stalder et. al. (2017)
#' set.seed(2018)
#' dat = simulateCC(ncase=500, ncontrol=500, beta0=-4.165,
#'                  betaG_SNP=c(log(1.2), log(1.2), 0, log(1.2), 0),
#'                  betaE_bin=log(1.5),
#'                  betaGE_SNP_bin=c(log(1.3), 0, 0, log(1.3), 0),
#'                  MAF=c(0.1, 0.3, 0.3, 0.3, 0.1),
#'                  SNP_cor=0.7, E_bin_freq=0.5)
#'
#' # SPMLE with known population disease rate of 0.03 and asymptotic SE estimates
#' spmleCombo(D=D, G=G, E=E, pi1=0.03, data=dat, nboot=0)
#' @export
spmleCombo = function(D, G, E, pi1, data, nboot=50, ncores=1, control=list(), startvals){
  if(nboot==0) {warning("nboot=0, using asymptotic standard error estimate, which has poor coverage properties")}

  ## Store the function call
  cl = match.call()

  ## Get argument names for D, G, and E
  Dname = substitute(D)
  Gname = substitute(G)
  Ename = substitute(E)

  ## If pi1 (partially) matches "rare", set pi1=0
  if(pmatch(x=substr(x=tolower(pi1), start=1, stop=4), table="rare", nomatch=FALSE)) {pi1=0}
  stopifnot(pi1>=0, pi1<1)  # throw an error if pi1 not in [0,1)

  ## Store the formula with user-provided variable names.  For consistency with estimators that accept formulas (like lm),
  ## set the formula environment to the calling environment (as if formula had been an argument).
  formula = formula(paste(as.character(as.expression(Dname)), "~", as.character(as.expression(Gname)), "*", as.character(as.expression(Ename))))
  attr(formula, ".Environment") = parent.frame()

  ## If no data.frame was supplied, set data to the calling environment
  if(missing(data)) {
    data = environment(formula)
  } else {                          # evaluate D, G, and E in data, if appropriate
    if(class(data)=="matrix") {data = as.data.frame(data)}
    D = with(data=data, eval(Dname))
    G = with(data=data, eval(Gname))
    E = with(data=data, eval(Ename))
  }

  ## Expand factor variables into dummies, if present
  if(class(G)=="factor") {G=model.matrix(~G)[,-1]}
  if(class(E)=="factor") {E=model.matrix(~E)[,-1]}

  ## Create model.frame (response, no interactions) and model.matrix (all predictors, including interactions)
  model_frame = model.frame(formula=formula, data=data)
  model_matrix = model.matrix(formula, model_frame)

  ## Fit the model with logistic regression
  logistic_fit = glm(formula, family=binomial(link='logit'), data=data)

  ## If starting values weren't provided, use logistic regression estimates
  if(missing(startvals)) {
    Omega_start = coef(logistic_fit)
  } else {
    Omega_start = startvals
  }

  ## If user-provided startval lacked names, add them
  if(is.null(names(Omega_start))) {names(Omega_start) = colnames(model_matrix)}

  ## Set control parameters
  con = list(trace=0, use_hess=TRUE, max_grad_tol=0.001, num_retries=2)
  con[(names(control))] = control

  ## Remove missing values and convert G and E to matrices
  complete = complete.cases(D, G, E)
  D = D[complete]
  G = as.matrix(as.matrix(G)[complete,])
  E = as.matrix(as.matrix(E)[complete,])

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  G = as.matrix(G)
  E = as.matrix(E)

  ## Model degrees of freedom
  df_model = length(Omega_start)

  ##### Treat G & E as in Stalder et al. 2017
  spmle_E = spmle(D=D, G=G, E=E, pi1=pi1, data=data, control=con, swap=FALSE, startvals=Omega_start)
  ##### Swap G & E
  spmle_G = spmle(D=D, G=G, E=E, pi1=pi1, data=data, control=con, swap=TRUE, startvals=Omega_start)

  ## Define matrices
  Omega_all = c(spmle_E$coefficients, spmle_G$coefficients)
  X_mat = rbind(diag(df_model), diag(df_model))
  zero_mat = matrix(0, nrow=df_model, ncol=df_model)

  ## Block diagonal with both inverse hessians
  H_all = -rbind(cbind(spmle_E$H_inv,zero_mat), cbind(zero_mat,spmle_G$H_inv))

  ## Sigma block matrix
  Sigma_EE  = spmle_E$Sigma
  Sigma_GG  = spmle_G$Sigma
  Sigma_EG  = ((ncontrol-1)*cov(spmle_E$zeta0, spmle_G$zeta0) + (ncase-1)*cov(spmle_E$zeta1, spmle_G$zeta1))/n
  Sigma_GE  = t(Sigma_EG)
  Sigma_all = rbind(cbind(Sigma_EE,Sigma_EG),cbind(Sigma_GE,Sigma_GG))

  ## Lambda (asymptotic covariance) matrix
  Lambda_all = H_all %*% Sigma_all %*% t(H_all)
  Lambda_all_inv = solve(Lambda_all)

  ## Asymptotic covariance matrix of optimum combination
  Lambda_combo = solve(t(X_mat) %*% Lambda_all_inv %*% X_mat)
  combo_asy_SE = sqrt(diag(Lambda_combo)/n)
  names(combo_asy_SE) = names(Omega_start)

  ## Optimum combination
  combo_par = as.vector(Lambda_combo %*% t(X_mat) %*% Lambda_all_inv %*% Omega_all)
  names(combo_par) = names(Omega_start)

  ## Calculate predictions
  linear_predictors = model_matrix %*% combo_par      # logistic scale
  fitted_values = plogis(q=linear_predictors)         # probability scale

  ## Pearson residuals, total deviance, and df resid & null
  pearson_resid = (D-fitted_values) / sqrt(fitted_values*(1-fitted_values))
  df_resid = n - df_model
  null_deviance = sum(-2*log(abs(1-D-mean(D))))
  df_null = n - 1

  ## Compile asymptotic results into a list.  Use glm-object naming conventions.
  combo_est = list(coefficients      = combo_par,
                   pi1               = pi1,
                   SE                = combo_asy_SE,
                   cov               = Lambda_combo/n,
                   spmle_E           = spmle_E,
                   spmle_G           = spmle_G,
                   glm_fit           = logistic_fit,
                   call              = cl,
                   formula           = formula,
                   data              = data,
                   model             = model_frame,
                   terms             = attr(model_frame, "terms"),
                   linear.predictors = linear_predictors,
                   fitted.values     = fitted_values,
                   residuals         = pearson_resid,
                   null.deviance     = null_deviance,
                   df.residual       = df_resid,
                   df.null           = df_null,
                   rank              = df_model,
                   nobs              = n,
                   ncase             = ncase,
                   ncontrol          = ncontrol)

  ## Set class spmle_combo
  class(combo_est) = c("spmle")

  ## Calculate bootstrap SE if nboot > 0
  if(nboot > 0) {
    if(con$trace > 0) {cat("Begin bootstrap\n")}

    ## Split the data into controls & cases for balanced resampling (same # of controls & cases in bootstraps as in original data)
    dat_control = list(D=D[D==0], G=G[D==0,], E=E[D==0,])
    dat_case = list(D=D[D==1], G=G[D==1,], E=E[D==1,])

    ## Within cases & controls, we use a balanced bootstrap: good asymptotic properties, easy to parallelize
    ## First create long lists of indices.  To create bootstrapped samples, read off this list of indices
    index_control = sample(rep(1:ncontrol, nboot))
    index_case = sample(rep(1:ncase, nboot))

    ## Parallelize bootstrap is ncores > 1
    if(ncores > 1) {
      ## Keep the user from specifying too many cores
      ncores = min(ncores, parallel::detectCores(logical=FALSE))

      ## Non-Windows operating systems can use mclapply with forking (faster)
      if(Sys.info()[1]!="Windows"){
        boot_ests = do.call(rbind, parallel::mclapply(X=seq_len(nboot), FUN=combo_boot, mc.cores=ncores,
                                            index_control=index_control, index_case=index_case,
                                            dat_control=dat_control, dat_case=dat_case, pi1=pi1,
                                            ncontrol=ncontrol, ncase=ncase, con=con))
      } else {  # windows uses PSOCK clusters
        cl = parallel::makeCluster(ncores)
        boot_ests = t(parallel::parSapply(cl=cl, X=seq_len(nboot), FUN=combo_boot,
                                          index_control=index_control, index_case=index_case,
                                          dat_control=dat_control, dat_case=dat_case, pi1=pi1,
                                          ncontrol=ncontrol, ncase=ncase, con=con))
        parallel::stopCluster(cl)
      }
    } else {  # otherwise, compute serially
      boot_ests = t(sapply(X=seq_len(nboot), FUN=combo_boot,
                           index_control=index_control, index_case=index_case,
                           dat_control=dat_control, dat_case=dat_case, pi1=pi1,
                           ncontrol=ncontrol, ncase=ncase, con=con))
    }  # end parallel bootstrap if()

    ## Calculate bootstrap covariance matrix and SE
    boot_cov = cov(boot_ests)
    boot_SE = sqrt(diag(boot_cov))

    ## Replace asymptotic covariance estimate with bootstrap
    combo_est$SE         = boot_SE
    combo_est$cov        = boot_cov
    combo_est$bootstraps = boot_ests
  }  # end bootstrap if()

  ## Return fitted model
  return(combo_est)
}


## Internal bootstrap function - can be called in parallel
combo_boot = function(b, index_control, index_case, dat_control, dat_case, pi1, ncontrol, ncase, con) {
  if(con$trace>0) {cat("Bootstrap", b, "\n")}

  ## Create bootstraped data set
  D_boot = c(rep(0,ncontrol),rep(1,ncase))
  G_boot = rbind(as.matrix(as.matrix(dat_control$G)[index_control[((b-1)*ncontrol+1):(b*ncontrol)],]), as.matrix(as.matrix(dat_case$G)[index_case[((b-1)*ncase+1):(b*ncase)],]))
  E_boot = rbind(as.matrix(as.matrix(dat_control$E)[index_control[((b-1)*ncontrol+1):(b*ncontrol)],]), as.matrix(as.matrix(dat_case$E)[index_case[((b-1)*ncase+1):(b*ncase)],]))

  ## Calculate both normal & swapped SPMLE with asymptotic SE
  spmle_E = spmle(D=D_boot, G=G_boot, E=E_boot, pi1=pi1, swap=FALSE, control=con)
  spmle_G = spmle(D=D_boot, G=G_boot, E=E_boot, pi1=pi1, swap=TRUE, control=con)

  ## Define matrices
  n = ncase + ncontrol
  df_model = length(spmle_E$coefficients)
  Omega_all = c(spmle_E$coefficients, spmle_G$coefficients)
  X_mat = rbind(diag(df_model), diag(df_model))
  zero_mat = matrix(0, nrow=df_model, ncol=df_model)

  ## Block diagonal with both inverse hessians
  H_all = -rbind(cbind(spmle_E$H_inv,zero_mat), cbind(zero_mat,spmle_G$H_inv))

  ## Sigma block matrix
  Sigma_EE  = spmle_E$Sigma
  Sigma_GG  = spmle_G$Sigma
  Sigma_EG  = ((ncontrol-1)*cov(spmle_E$zeta0, spmle_G$zeta0) + (ncase-1)*cov(spmle_E$zeta1, spmle_G$zeta1))/n
  Sigma_GE  = t(Sigma_EG)
  Sigma_all = rbind(cbind(Sigma_EE,Sigma_EG),cbind(Sigma_GE,Sigma_GG))

  ## Lambda (asymptotic covariance) matrix
  Lambda_all = H_all %*% Sigma_all %*% t(H_all)
  Lambda_all_inv = solve(Lambda_all)

  ## Asymptotic covariance matrix of optimum combination
  Lambda_combo = solve(t(X_mat) %*% Lambda_all_inv %*% X_mat)

  ## Optimum combination
  combo_par = as.vector(Lambda_combo %*% t(X_mat) %*% Lambda_all_inv %*% Omega_all)

  return(combo_par)
}


