


################################################################################ SPMLE

## Wrapper to calculate negative loglikelihood and gradient.  Report likelihood, store gradient
lik_fn = function(Omega, D, G, E, pi1) {
  lik_grad = neglikgrad(Omega=Omega, D=D, G=G, E=E, pi1=pi1) # neglikgrad is c++ code that calculates the negative loglikelihood and gradient
  assign(x="grad_val", value=lik_grad[["gradient"]], pos=parent.frame())  # Store the gradient in the calling environment to retrieve it when ucminf tries to calculate the gradient
  # grad_val <<- lik_grad[[2]]  # Store the gradient in .GlobalEnv to retrieve it when ucminf tries to calculate the gradient
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

  ## Check for convergance and rerun optimization if necessary
  if(!is.finite(spmle_max$info[1]) | spmle_max$info[1]>con$max_grad_tol) {
    ## Scale starting values by 1/SD
    scale_vec = 1/c(1, apply(model.matrix(~0+G*E), 2, sd))

    ## Loop through different starting values
    for(j in seq_len(con$num_retries)) {
      if(con$trace>-1) {cat("ucminf retry", j, "of", con$num_retries, "\n")}

      ## If failure happened when preconditioning with the hessian, try without (same start values the first time)
      if(!is.null(ucminf_con$invhessian.lt)) {
        ucminf_con$invhessian.lt = NULL
        spmle_max = ucminf::ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)
      } else {  # otherwise use scaled normal startvals (preserve names from Omega_start)
        spmle_max = ucminf::ucminf(par=setNames(rnorm(n=length(Omega_start), sd=scale_vec), names(Omega_start)), fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)
      }

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
#' @param D a binary vector of disease status (1=case, 0=control).
#' @param G a vector or matrix (if multivariate) containing genetic data. Can be continuous, discrete, or a combination.
#' @param E a vector or matrix (if multivariate) containing environmental data. Can be continuous, discrete, or a combination.
#' @param pi1 the population disease rate, a scalar in [0, 1).  Using \code{pi1=0} is the rare disease approximation.
#' @param data an optional data frame, list, or environment (or object coercible by \code{\link[base]{as.data.frame}}
#'   to a data frame) containing the variables in the model.  If not found in data, the variables are taken from
#'   the environment from which \code{spmle} is called.
#' @param control a list of control parameters that allow the user to control the optimization algorithm.  See 'Details'.
#' @param swap a logical scalar rarely of interest to the end user.  Dependence on the distributions of G and E
#'   are removed using different methods; this switch swaps them to produce a symmetric estimator with identical
#'   properties to the SPMLE.  Default \code{FALSE}.
#' @param startvals an optional numeric vector of coefficient starting values for optimization.  Usually left blank,
#'   in which case logistic regression estimates are used as starting values.
#' @return an object of class \code{spmle}.
#' @seealso \code{\link{simulate_complex}} to simulate data
#' @examples
#' # Simulation from Table 1 in Stalder et. al. (2017)
#' dat = simulate_complex(ncase=1000,
#'                        ncontrol=1000,
#'                        beta0=-4.14,
#'                        betaG_SNP=c(log(1.2), log(1.2), 0, log(1.2), 0),
#'                        betaE_bin=log(1.5),
#'                        betaGE_SNP_bin=c(log(1.3), 0, 0, log(1.3), 0),
#'                        MAF=c(0.1, 0.3, 0.3, 0.3, 0.1),
#'                        SNP_cor=0.7,
#'                        E_bin_freq=0.5)
#'
#' # SPMLE with known population disease rate of 0.03
#' spmle(D=dat$D, G=dat$G, E=dat$E, pi1=0.03)
#'
#' # Simulation with a single SNP and a single binary environmental variable.
#' # True population disease rate in this simulation is 0.03.
#' # This simulation scenario was used in the Supplementary Material of Stalder et. al. (2017)
#' # to compare performance against the less flexible method of Chatterjee and Carroll (2005),
#' # which is available as the function as snp.logistic in the Bioconductor package CGEN.
#' dat = simulate_complex(ncase=100,
#'                        ncontrol=100,
#'                        beta0=-3.77,
#'                        betaG_SNP=log(1.2),
#'                        betaE_bin=log(1.5),
#'                        betaGE_SNP_bin=log(1.3),
#'                        MAF=0.1)
#'
#' # SPMLE using the rare disease assumption
#' #and with bootstrap SE, no tracing, and no hessian preconditioning.
#' spmle(D=dat$D, G=dat$G, E=dat$E, pi1=0, control=list(nboot=100, trace=0, use_hess=TRUE))
#' @export
spmle = function(D, G, E, pi1, data, control=list(), swap=FALSE, startvals){
  ## Store the function call
  cl = match.call()

  ## Get argument names for D, G, and E
  Dname = substitute(D)
  Gname = substitute(G)
  Ename = substitute(E)

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

  ## If starting values weren't provided, use logistic regression estimates
  if(missing(startvals)) {
    Omega_start = coef(glm(formula, family=binomial(link='logit'), data=data))
  } else {
    Omega_start = startvals
  }

  ## If user-provided startval lacked names, add them
  if(is.null(names(Omega_start))) {names(Omega_start) = colnames(model_matrix)}

  ## Set control parameters
  con = list(trace=0, use_hess=TRUE, max_grad_tol=0.001, num_retries=2)
  con[(names(control))] = control

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  G = as.matrix(G)
  E = as.matrix(E)

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
  if(con$use_hess==TRUE) {con$invhessian.lt = solve(neghess(Omega=Omega_start, D=D, G=G, E=E, pi1=pi1))[lower.tri(diag(length(Omega_start)),diag=TRUE)]}

  ## Calculate SPMLE for a given disease rate
  spmle_max = maximize_spmle(Omega_start=Omega_start[swap_order], D=D, G=G, E=E, pi1=pi1, control=con)

  ## Asymptotic SE for known disease rate
  hess_zeta = hesszeta(Omega=spmle_max$par, D=D, G=G, E=E, pi1=pi1)
  Sigma = ((ncontrol-1)*cov(hess_zeta$zeta0) + (ncase-1)*cov(hess_zeta$zeta1))/n  # (ncontrol-1) to correct denominator because cov uses the "sample covariance matrix" denominator of n-1
  H_inv = -n*solve(hess_zeta$hessian)    # = (-hessian/n)^-1 = (Gamma_1 - Gamma_2)^-1
  Lambda = H_inv %*% Sigma %*% t(H_inv)  # covar matrix of sqrt(n) * OmegaHat
  SE_asy = sqrt(diag(Lambda)/n)

  ## Calculate predictions
  linear_predictors = model_matrix %*% spmle_max$par  # logistic scale
  fitted_values = plogis(q=linear_predictors)         # probability scale

  ## Deviance residuals and total deviance
  deviance_resid = -2*log(abs(1-D-fitted_values))
  total_deviance = sum(deviance_resid)

  ## Compile results into a list.  Use glm-object naming conventions.  If we swapped G & E, change back now
  spmle_est = list(coefficients = spmle_max$par[reverse_swap_order],
                   SE    = SE_asy[reverse_swap_order],
                   cov   = Lambda[reverse_swap_order,reverse_swap_order]/n,
                   H_inv = H_inv[reverse_swap_order,reverse_swap_order],
                   Sigma = Sigma[reverse_swap_order,reverse_swap_order],
                   zeta0 = hess_zeta$zeta0[, reverse_swap_order],
                   zeta1 = hess_zeta$zeta1[, reverse_swap_order],
                   ucminf = spmle_max,
                   call  = cl,
                   formula = formula,
                   data = data,
                   model = model_frame,
                   linear.predictors = linear_predictors,
                   fitted.values = fitted_values,
                   residuals = deviance_resid,
                   deviance = total_deviance,
                   aic = NA)

  ## Return an object of class spmle, which inherits from glm and lm
  class(spmle_est) = c("spmle")#, "glm", "lm")
  return(spmle_est)
}


################################################################################ Combo SPMLE

#' Calculate seven estimators with asymptotic SEs: logistic, SPMLE, swap, Optimal Combo
#'
#' @param D a binary vector of disease status (1=case, 0=control).
#' @param G a vector or matrix (if multivariate) containing genetic data. Can be continuous, discrete, or a combination.
#' @param E a vector or matrix (if multivariate) containing environmental data. Can be continuous, discrete, or a combination.
#' @param pi1 the population disease rate, a scalar in [0, 1).  Using \code{pi1=0} is the rare disease approximation.
#' @param control a list of control parameters.
#'     \code{nboot} number of bootstraps.  Default \code{0}.
#'     \code{trace} is a scalar.  If >-1, tracing information is produced.  Default \code{0}.
#'     \code{use_hess} logical: precondition optimization with the hessian.  Default \code{FALSE}.
#' @return A list with a matrix of estimates and SEs and a bunch of other stuff
#' @seealso \code{\link{simulate_complex}} to simulate data
#' @export
combo_asymp = function(D, G, E, pi1, control=list()){
  ## Set control parameters
  con = list(nboot=0, trace=0, use_hess=TRUE)
  con[(names(control))] = control

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  G = as.matrix(G)
  E = as.matrix(E)

  ## Use Logistic estimates as starting values
  logistic_fit = glm(D~G*E, family=binomial(link='logit'))
  logistic_est = t(coef(summary(logistic_fit))[,c(1,2)])
  Omega_start  = logistic_est[1,]
  length_Omega = length(Omega_start)

  ##### Treat G & E as in Stalder et al. 2017
  spmle_G_asy = spmle(startvals=Omega_start, D=D, G=G, E=E, pi1=pi1, swap=FALSE, control=con)
  ##### Swap G & E
  spmle_E_asy = spmle(startvals=Omega_start, D=D, G=G, E=E, pi1=pi1, swap=TRUE, control=con)

  ## Define matrices
  Omega_all = c(spmle_G_asy$coefficients, spmle_E_asy$coefficients)
  X_mat = rbind(diag(length_Omega), diag(length_Omega))
  zero_mat = matrix(0, nrow=length_Omega, ncol=length_Omega)

  ## Block diagonal with both inverse hessians
  H_all = -rbind(cbind(spmle_G_asy$H_inv,zero_mat), cbind(zero_mat,spmle_E_asy$H_inv))

  ## Sigma block matrix
  Sigma_GG  = spmle_G_asy$Sigma
  Sigma_EE  = spmle_E_asy$Sigma
  Sigma_GE  = ((ncontrol-1)*cov(spmle_G_asy$zeta0, spmle_E_asy$zeta0) + (ncase-1)*cov(spmle_G_asy$zeta1, spmle_E_asy$zeta1))/n
  Sigma_EG  = t(Sigma_GE)
  Sigma_all = rbind(cbind(Sigma_GG,Sigma_GE),cbind(Sigma_EG,Sigma_EE))

  ## Lambda (asymptotic covariance) matrix
  Lambda_all = H_all %*% Sigma_all %*% t(H_all)
  Lambda_all_inv = solve(Lambda_all)

  ## Covariance matrix of optimum combination
  Lambda_combo = solve(t(X_mat) %*% Lambda_all_inv %*% X_mat)
  combo_SE = sqrt(diag(Lambda_combo)/n)

  ## Optimum combination
  combo_par = as.vector(Lambda_combo %*% t(X_mat) %*% Lambda_all_inv %*% Omega_all)

  ## All seven estimates
  asy_ests = rbind(logistic_est,
                   spmle_G_asy$coefficients, spmle_G_asy$SE,
                   spmle_E_asy$coefficients, spmle_E_asy$SE,
                   combo_par, combo_SE
                   )
  rownames(asy_ests) = c("logistic_par", "logistic_SE",
                         "spmle_G_par", "spmle_G_SE",
                         "spmle_E_par", "spmle_E_SE",
                         "combo_par", "combo_SE"
                         )

  ## Return estimates
  return(list(ests              = asy_ests,
  						Lambda_all        = Lambda_all,
  						Sigma_GG          = Sigma_GG,
  						Sigma_EE          = Sigma_EE,
  						Sigma_GE          = Sigma_GE,
  						Sigma_all         = Sigma_all,
  						Lambda_combo      = Lambda_combo,
  						spmle_G_asy_zeta0 = spmle_G_asy$zeta0,
  						spmle_G_asy_zeta1 = spmle_G_asy$zeta1,
  						spmle_E_asy_zeta0 = spmle_E_asy$zeta0,
  						spmle_E_asy_zeta1 = spmle_E_asy$zeta1,
  						spmle_G_asy_H     = spmle_G_asy$H,
  						spmle_E_asy_H     = spmle_E_asy$H
  ))
}

#' \code{nboot} is the number of bootstrap samples to be used when calculating the bootstrap SE.  Setting \code{nboot=0},
#'  the default, will skip the bootstrap and calculate only asymptotic SE.
#' Function to calculate bootstrapped SE for spmle_G & spmle_E estimates, as well as bootstrapped Combo estimate
#'
#' @param D a binary vector of disease status (1=case, 0=control).
#' @param G a vector or matrix (if multivariate) containing genetic data. Can be continuous, discrete, or a combination.
#' @param E a vector or matrix (if multivariate) containing environmental data. Can be continuous, discrete, or a combination.
#' @param pi1 the population disease rate, a scalar in [0, 1).  Using \code{pi1=0} is the rare disease approximation.
#' @param control a list of control parameters.
#'     \code{nboot} number of bootstraps.  Default \code{0}.
#'     \code{trace} is a scalar.  If >-1, tracing information is produced.  Default \code{0}.
#'     \code{use_hess} logical: precondition optimization with the hessian.  Default \code{FALSE}.
#' @param spmle_G_par,spmle_E_par Parameter estimates generated from \code{\link{combo_asymp}}
#' @return A list with a matrix of estimates and SEs and a bunch of other stuff
#' @seealso \code{\link{simulate_complex}} to simulate data
#' @export
combo_boot = function(D, G, E, pi1, spmle_G_par, spmle_E_par, control=list()) {
  ## Set control parameters
  con = list(nboot=0, trace=0, use_hess=TRUE)
  con[(names(control))] = control

  ## Correct too-small bootstrap sample
  if(con$nboot!=0 && con$nboot<2) {
    warning("To use the bootstrap, nboot must be >= 2.  Skipping bootstrap.")
    con$nboot = 0
  }
  nboot = con$nboot

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  length_Omega = length(spmle_G_par)
  G = as.matrix(G)
  E = as.matrix(E)
  nG = NCOL(G)
  nE = NCOL(E)
  reverse_swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))

  ############### Bootstrap ####################
  ## Split the data into controls & cases for balanced resampling (same # of controls & cases in bootstraps as in original data)
  dat_control = list(D=D[D==0], G=G[D==0,], E=E[D==0,])
  dat_case = list(D=D[D==1], G=G[D==1,], E=E[D==1,])

  ## Within cases & controls, we use a balanced bootstrap: good asymptotic properties, easy to parallelize
  ## First create long lists of indices.  To create bootstrapped samples, read off this list of indices
  index_control = sample(rep(1:ncontrol, nboot))
  index_case = sample(rep(1:ncase, nboot))

  spmle_G_boot = spmle_E_boot = symple_fusion_boot = matrix(nrow=nboot, ncol=length_Omega)  # Set empty matrix to store bootstrap results
  for(b in 1:nboot) { # bootstrap loop
    if((b %% ceiling(sqrt(nboot))) == 0 && con$trace > -1) {
      updatetxt = paste0("Bootstrap ", b,"\n")
      cat(updatetxt)
    }

    ## Create bootstraped data set
    D_boot = c(rep(0,ncontrol),rep(1,ncase))
    G_boot = rbind(as.matrix(as.matrix(dat_control$G)[index_control[((b-1)*ncontrol+1):(b*ncontrol)],]), as.matrix(as.matrix(dat_case$G)[index_case[((b-1)*ncase+1):(b*ncase)],]))
    E_boot = rbind(as.matrix(as.matrix(dat_control$E)[index_control[((b-1)*ncontrol+1):(b*ncontrol)],]), as.matrix(as.matrix(dat_case$E)[index_case[((b-1)*ncase+1):(b*ncase)],]))

    ## Logistic regression
    logistic_fit = glm(D_boot~G_boot*E_boot, family=binomial(link='logit'))

    ## Calculate both normal & swapped SPMLE with asymptotic SE
   	spmle_G_asymp_boot = spmle(startvals=coef(logistic_fit), D=D_boot, G=G_boot, E=E_boot, pi1=pi1, swap=FALSE, control=con)
   	spmle_E_asymp_boot = spmle(startvals=coef(logistic_fit), D=D_boot, G=G_boot, E=E_boot, pi1=pi1, swap=TRUE, control=con)
   	spmle_G_boot[b,] = spmle_G_asymp_boot$coefficients
   	spmle_E_boot[b,] = spmle_E_asymp_boot$coefficients

    ### Fusion (combo) asymptotic (to calculate bootstrap SE even when using asymptotic point estimate)
    ## Define matrices
    Omega_all_boot = c(spmle_G_asymp_boot$coefficients, spmle_E_asymp_boot$coefficients)
    X_mat = rbind(diag(length_Omega), diag(length_Omega))
    zero_mat = matrix(0, nrow=length_Omega, ncol=length_Omega)

    ## Block diagonal with both inverse hessians
    H_all_boot = -rbind(cbind(spmle_G_asymp_boot$H_inv, zero_mat), cbind(zero_mat, spmle_E_asymp_boot$H_inv))

    ## Sigma block matrix
    Sigma_GG_boot  = spmle_G_asymp_boot$Sigma
    Sigma_EE_boot  = spmle_E_asymp_boot$Sigma
    Sigma_GE_boot  = ((ncontrol-1)*cov(spmle_G_asymp_boot$zeta0, spmle_E_asymp_boot$zeta0) + (ncase-1)*cov(spmle_G_asymp_boot$zeta1, spmle_E_asymp_boot$zeta1))/n
    Sigma_EG_boot  = t(Sigma_GE_boot)
    Sigma_all_boot = rbind(cbind(Sigma_GG_boot, Sigma_GE_boot),cbind(Sigma_GE_boot, Sigma_EE_boot))

    ## Lambda (asymptotic covariance) matrix
    Lambda_all_boot = H_all_boot %*% Sigma_all_boot %*% t(H_all_boot)
    Lambda_all_inv_boot = solve(Lambda_all_boot)

    ## Optimum combination
    Lambda_combo_boot = solve(t(X_mat) %*% Lambda_all_inv_boot %*% X_mat)
    symple_fusion_boot[b,] = as.vector(Lambda_combo_boot %*% t(X_mat) %*% Lambda_all_inv_boot %*% Omega_all_boot)

  }  ## End bootstrap loop

  ### Calculate bootstrap SE
  symple_fusion_asymp_boot_SE = apply(symple_fusion_boot, 2, FUN=function(x) sd(x, na.rm = TRUE))

  ### Calculate bootstrap combo estimator
  Omega_all = c(spmle_G_par, spmle_E_par)
  X_mat = rbind(diag(length_Omega), diag(length_Omega))

  ## Return estimates and Lambda matrix
  return(list(boot_SE = symple_fusion_asymp_boot_SE,
  						spmle_G_boot = spmle_G_boot,
  						spmle_E_boot = spmle_E_boot,
  						symple_fusion_boot = symple_fusion_boot
  						))
}


################################################################################ Composite SPMLE

## Sums the neg loglik from SPMLE and swapped likelihoods
composite_lik = function(Omega, D, G, E, pi1) {
  nG = NCOL(G)
  nE = NCOL(E)
  reverse_swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
  swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))
  spmle_G_lik = neglikgrad(Omega=Omega, D=D, G=G, E=E, pi1=pi1)  # neg loglik from spmle
  spmle_E_lik = neglikgrad(Omega=Omega[swap_order], D=D, G=E, E=G, pi1=pi1)  # neg loglik from swap
  assign(x="grad_val", value=spmle_G_lik$gradient + spmle_E_lik$gradient[reverse_swap_order], pos=parent.frame())  # Store the gradient in the calling environment to retrieve it when ucminf tries to calculate the gradient
  # grad_val <<- spmle_G_lik$gradient + spmle_E_lik$gradient[reverse_swap_order] # Store the gradient in .GlobalEnv to retrieve it when ucminf tries to calculate the gradient
  return(spmle_G_lik$objective + spmle_E_lik$objective)
}


## Calculate the composite SPMLE (maximize sum of normal SPMLE and swapped likelihoods)
composite_est = function(Omega_start, D, G, E, pi1, control=list()) {
  ## Set control parameters
  con = list(nboot=0, trace=0, use_hess=TRUE, max_grad_tol=0.001, num_retries=21)
  con[(names(control))] = control

  nG = NCOL(G)
  nE = NCOL(E)
  reverse_swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
  swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))

  if(con$use_hess==TRUE) {
    ## Calculate hessian
    hg = neghess(Omega=Omega_start, D=D, G=G, E=E, pi1=pi1)
    he = neghess(Omega=Omega_start[swap_order], D=D, G=E, E=G, pi1=pi1)
    H0 = (hg + he[reverse_swap_order,reverse_swap_order])

    ## Lower triangle of inverse hessian (use identity matrix if hessian is not PD)
    invhessian.lt = tryCatch(solve(H0), error=function(x) diag(length(Omega_start)))[lower.tri(diag(length(Omega_start)),diag=TRUE)]

    ## Maximize composite likelihood
    composite_spmle = ucminf::ucminf(par=Omega_start, fn=composite_lik, gr=grad_fn, control=list(trace=con$trace, invhessian.lt=invhessian.lt), D=D, G=G, E=E, pi1=pi1)
  } else {
    composite_spmle = ucminf::ucminf(par=Omega_start, fn=composite_lik, gr=grad_fn, control=list(trace=con$trace), D=D, G=G, E=E, pi1=pi1)
  }

  ## Check for convergance, and rerun optimization if necessary
  if(!is.finite(composite_spmle$info[1]) | composite_spmle$info[1]>con$max_grad_tol) {
    startvals=cbind(rep(0,length(Omega_start)), # create a matrix of starting values
                    matrix(runif(length(Omega_start)*ceiling((con$num_retries-1)/4), min=-1, max=1), nrow=length(Omega_start)) + Omega_start,
                    matrix(rnorm(length(Omega_start)*ceiling((con$num_retries-1)/4)), nrow=length(Omega_start)),
                    matrix(rt(length(Omega_start)*ceiling((con$num_retries-1)/2), df=2), nrow=length(Omega_start)))
    for(j in 1:ncol(startvals)) { # retry optimization
      updatetxt = paste0("ucminf retry ",j,"\n")
      if(con$trace>-1) {cat(updatetxt)}
      composite_spmle = ucminf::ucminf(par=startvals[,j], fn=composite_lik, gr=grad_fn, control=list(trace=con$trace), D=D, G=G, E=E, pi1=pi1)
      if(is.finite(composite_spmle$info[1]) & composite_spmle$info[1]<con$max_grad_tol) break  # stop once we have convergence
      if(j == ncol(startvals)) stop("ucminf failed to converge")
    }  # end of j loop
  }
  return(composite_spmle)
}


## Calculate asymptotic SE for NPMLE
composite_asymp = function(D, G, E, pi1, Omega_start=NULL, control=list()){
  ## Set control parameters
  con = list(nboot=0, trace=0, use_hess=TRUE)
  con[(names(control))] = control

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  G = as.matrix(G)
  E = as.matrix(E)
  nG = NCOL(G)
  nE = NCOL(E)

  reverse_swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
  swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))

  ## Use Logistic estimates as starting values if they weren't provided
  if(is.null(Omega_start)){
    logistic_fit = glm(D~G*E, family=binomial(link='logit'))
    logistic_est = t(coef(summary(logistic_fit))[,c(1,2)])
    Omega_start  = logistic_est[1,]
  }

  ## Calculate SPMLE for a given disease rate
  composite_par = composite_est(Omega_start=Omega_start, D=D, G=G, E=E, pi1=pi1, control=con)$par

  ### Asymptotic SE for known disease rate
  ## First calculate hessians (Gamma 1 & 2), and zetas for SPMLE_G and SPMLE_E
  hess_zeta_G = hesszeta(Omega=composite_par, D=D, G=G, E=E, pi1=pi1)
  hess_zeta_E = hesszeta(Omega=composite_par[swap_order], D=D, G=E, E=G, pi1=pi1)

  ## Combined hessian
  H = (hess_zeta_G$hessian + hess_zeta_E$hessian[reverse_swap_order,reverse_swap_order])/2
  H_inv = -n*solve(H)             # = -(hessian/n)^(-1) = (Gamma_1 - Gamma_2)^-1

  ## Combined zetas
  Zeta0 = (hess_zeta_G$zeta0 + hess_zeta_E$zeta0[, reverse_swap_order])/2
  Zeta1 = (hess_zeta_G$zeta1 + hess_zeta_E$zeta1[, reverse_swap_order])/2
  Sigma = ((ncontrol-1)*cov(Zeta0) + (ncase-1)*cov(Zeta1))/n

  ## Asymptotic SE for the composite estimator follows the same pattern as the standard SPMLE
  Lambda = H_inv %*% Sigma %*% t(H_inv)  # covar matrix of sqrt(n) * OmegaHat
  SE_asy = sqrt(diag(Lambda)/n)

  ## Return results
  return(list(coefficients = composite_par,
              SE    = SE_asy,
              cov   = Lambda,
              H_inv = H_inv,
              Sigma = Sigma,
              zeta0 = Zeta0,
              zeta1 = Zeta1,
              H     = H)
  )
}
