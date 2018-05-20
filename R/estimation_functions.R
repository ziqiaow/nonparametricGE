


################################################################################ SPMLE

## Calculate negative loglikelihood and gradient.  Report likelihood, store gradient
lik_fn = function(Omega, D, G, E, pi1) {
  lik_grad = neglikgrad(Omega=Omega, D=D, G=G, E=E, pi1=pi1) # neglikgrad is c++ code that calculates the negative loglikelihood and gradient
  assign(x="grad_val", value=lik_grad[[2]], pos=parent.frame())  # Store the gradient in the calling environment to retrieve it when ucminf tries to calculate the gradient
  # grad_val <<- lik_grad[[2]]  # Store the gradient in .GlobalEnv to retrieve it when ucminf tries to calculate the gradient
  return(lik_grad[[1]])
}

## Report the gradient
grad_fn = function(Omega, D, G, E, pi1) return(get(x="grad_val", pos=parent.frame())) # return the previously stored gradient

## Calculate the SPMLE (optionally precondition with hessian)
SPMLE_fun = function(Omega_start, D, G, E, pi1, control=list()) {
  ## Set control parameters for SPMLE estimation
  con = list(trace=0, max_grad_tol=0.001, num_tries=21)
  con[(names(control))] = control

  ## Set control parameters for ucminf (using values from con where appropriate)
  ucminf_con = list(grtol=1e-06, xtol=1e-12, stepmax=1, maxeval=500,
                    grad="forward", gradstep=c(1e-06, 1e-08), invhessian.lt=NULL)
  ucminf_con[names(ucminf_con) %in% names(con)] = con[names(con) %in% names(ucminf_con)]

  ## Optimize with ucminf.  Here Omega_start is the starting values for optimization and invhessian.lt is the lower triangle of the inverse hessian
  SPMLE = ucminf::ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)

  ## When ucminf works, it works brilliantly.  But it has a nasty habit of declaring
  ## convergence before actually converging.  Here we check for convergance
  ## and rerun optimization if necessary.
  ##
  ## max_grad_tol = maximum allowable gradient
  ## num_tries = number of times to retry optimization
  if(!is.finite(SPMLE$info[1]) | SPMLE$info[1]>con$max_grad_tol) {
    ucminf_con$invhessian.lt = NULL  # don't use hessian during retries
    startvals = cbind(rep(0,length(Omega_start)), # create a matrix of starting values
                      matrix(runif(length(Omega_start)*ceiling((con$num_tries-1)/4), min=-1, max=1), nrow=length(Omega_start)) + Omega_start,
                      matrix(rnorm(length(Omega_start)*ceiling((con$num_tries-1)/4)), nrow=length(Omega_start)),
                      matrix(rt(length(Omega_start)*ceiling((con$num_tries-1)/2), df=2), nrow=length(Omega_start)))
    for(j in 1:ncol(startvals)) { # retry optimization
      updatetxt = paste0("UCMINF retry ",j,"\n")
      if(con$trace>-1) cat(updatetxt)
      SPMLE = ucminf::ucminf(par=startvals[,j], fn=lik_fn, gr=grad_fn, control=list(trace=con$trace), D=D, G=G, E=E, pi1=pi1)
      if(is.finite(SPMLE$info[1]) & SPMLE$info[1]<con$max_grad_tol) break  # stop once we have convergence
      if(j == ncol(startvals)) stop("UCMINF failed to converge")
    }  # end of j loop
  }
  return(SPMLE$par)
}

## Calculate asymptotic SE for SPMLE
SPMLE_asymp = function(D, G, E, pi1, control=list(), swap=FALSE, Omega_start=NULL){
  ## Set control parameters
  con = list(nboot=0, trace=0, usehess=FALSE)
  con[(names(control))] = control

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  G = as.matrix(G)
  E = as.matrix(E)

  ## Use Logistic estimates as starting values if they weren't provided
  if(is.null(Omega_start)) Omega_start = coef(glm(D~G*E, family=binomial(link='logit')))

  ## If we're swapping G & E, make the change now
  if(swap==TRUE) {
    nG = NCOL(G)
    nE = NCOL(E)
    swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
    reverse_swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))
    temp = G
    G = E
    E = temp
  } else {
    swap_order = reverse_swap_order = seq_along(Omega_start)
  }

  ## Calculate SPMLE for a given disease rate
  if(con$usehess==TRUE) {
    H0 = neghess(Omega=Omega_start[reverse_swap_order], D=D, G=G, E=E, pi1=pi1)
  } else {
    H0 = NULL
  }

  SPMLE = SPMLE_fun(Omega_start=Omega_start[reverse_swap_order], D=D, G=G, E=E, pi1=pi1, H0=H0, control=con)

  ## Asymptotic SE for known disease rate
  hess_zeta = hesszeta(Omega=SPMLE, D=D, G=G, E=E, pi1=pi1)
  Sigma = ((ncontrol-1)*cov(hess_zeta$zeta0) + (ncase-1)*cov(hess_zeta$zeta1))/n  # (ncontrol-1) to correct denominator because cov uses the "sample covariance matrix" denominator of n-1
  H_inv = -n*solve(hess_zeta$hessian)    # = (-hessian/n)^-1 = (Gamma_1 - Gamma_2)^-1
  Lambda = H_inv %*% Sigma %*% t(H_inv)  # covar matrix of sqrt(n) * OmegaHat
  SE_asy = sqrt(diag(Lambda)/n)

  ## If we're swapping G & E, change back now
  return(list(par   = SPMLE[swap_order],
              SE    = SE_asy[swap_order],
              cov   = Lambda[swap_order,swap_order],
              H_inv = H_inv[swap_order,swap_order],
              Sigma = Sigma[swap_order,swap_order],
              zeta0 = hess_zeta$zeta0[, swap_order],
              zeta1 = hess_zeta$zeta1[, swap_order],
  						H     = hess_zeta$hessian[swap_order,swap_order])
  )
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
#'     \code{usehess} logical: precondition optimization with the hessian.  Default \code{FALSE}.
#' @return A list with a matrix of estimates and SEs and a bunch of other stuff
#' @seealso \code{\link{simulate_complex}} to simulate data
#' @export
combo_asymp = function(D, G, E, pi1, control=list()){
  ## Set control parameters
  con = list(nboot=0, trace=0, usehess=FALSE)
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
  SPMLE_G_asy = SPMLE_asymp(Omega_start=Omega_start, D=D, G=G, E=E, pi1=pi1, swap=FALSE, control=con)
  ##### Swap G & E
  SPMLE_E_asy = SPMLE_asymp(Omega_start=Omega_start, D=D, G=G, E=E, pi1=pi1, swap=TRUE, control=con)

  ## Define matrices
  Omega_all = c(SPMLE_G_asy$par, SPMLE_E_asy$par)
  X_mat = rbind(diag(length_Omega), diag(length_Omega))
  zero_mat = matrix(0, nrow=length_Omega, ncol=length_Omega)

  ## Block diagonal with both inverse hessians
  H_all = -rbind(cbind(SPMLE_G_asy$H_inv,zero_mat), cbind(zero_mat,SPMLE_E_asy$H_inv))

  ## Sigma block matrix
  Sigma_GG  = SPMLE_G_asy$Sigma
  Sigma_EE  = SPMLE_E_asy$Sigma
  Sigma_GE  = ((ncontrol-1)*cov(SPMLE_G_asy$zeta0, SPMLE_E_asy$zeta0) + (ncase-1)*cov(SPMLE_G_asy$zeta1, SPMLE_E_asy$zeta1))/n
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
                   SPMLE_G_asy$par, SPMLE_G_asy$SE,
                   SPMLE_E_asy$par, SPMLE_E_asy$SE,
                   combo_par, combo_SE
                   )
  rownames(asy_ests) = c("logistic_par", "logistic_SE",
                         "SPMLE_G_par", "SPMLE_G_SE",
                         "SPMLE_E_par", "SPMLE_E_SE",
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
  						SPMLE_G_asy_zeta0 = SPMLE_G_asy$zeta0,
  						SPMLE_G_asy_zeta1 = SPMLE_G_asy$zeta1,
  						SPMLE_E_asy_zeta0 = SPMLE_E_asy$zeta0,
  						SPMLE_E_asy_zeta1 = SPMLE_E_asy$zeta1,
  						SPMLE_G_asy_H     = SPMLE_G_asy$H,
  						SPMLE_E_asy_H     = SPMLE_E_asy$H
  ))
}


#' Function to calculate bootstrapped SE for SPMLE_G & SPMLE_E estimates, as well as bootstrapped Combo estimate
#'
#' @param D a binary vector of disease status (1=case, 0=control).
#' @param G a vector or matrix (if multivariate) containing genetic data. Can be continuous, discrete, or a combination.
#' @param E a vector or matrix (if multivariate) containing environmental data. Can be continuous, discrete, or a combination.
#' @param pi1 the population disease rate, a scalar in [0, 1).  Using \code{pi1=0} is the rare disease approximation.
#' @param control a list of control parameters.
#'     \code{nboot} number of bootstraps.  Default \code{0}.
#'     \code{trace} is a scalar.  If >-1, tracing information is produced.  Default \code{0}.
#'     \code{usehess} logical: precondition optimization with the hessian.  Default \code{FALSE}.
#' @param SPMLE_G_par,SPMLE_E_par Parameter estimates generated from \code{\link{combo_asymp}}
#' @return A list with a matrix of estimates and SEs and a bunch of other stuff
#' @seealso \code{\link{simulate_complex}} to simulate data
#' @export
combo_boot = function(D, G, E, pi1, SPMLE_G_par, SPMLE_E_par, control=list()) {
  ## Set control parameters
  con = list(nboot=0, trace=0, usehess=FALSE)
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
  length_Omega = length(SPMLE_G_par)
  G = as.matrix(G)
  E = as.matrix(E)
  nG = NCOL(G)
  nE = NCOL(E)
  swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))

  ############### Bootstrap ####################
  ## Split the data into controls & cases for balanced resampling (same # of controls & cases in bootstraps as in original data)
  dat_control = list(D=D[D==0], G=G[D==0,], E=E[D==0,])
  dat_case = list(D=D[D==1], G=G[D==1,], E=E[D==1,])

  ## Within cases & controls, we use a balanced bootstrap: good asymptotic properties, easy to parallelize
  ## First create long lists of indices.  To create bootstrapped samples, read off this list of indices
  index_control = sample(rep(1:ncontrol, nboot))
  index_case = sample(rep(1:ncase, nboot))

  SPMLE_G_boot = SPMLE_E_boot = symple_fusion_boot = matrix(nrow=nboot, ncol=length_Omega)  # Set empty matrix to store bootstrap results
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
   	SPMLE_G_asymp_boot = SPMLE_asymp(Omega_start=coef(logistic_fit), D=D_boot, G=G_boot, E=E_boot, pi1=pi1, swap=FALSE, control=con)
   	SPMLE_E_asymp_boot = SPMLE_asymp(Omega_start=coef(logistic_fit), D=D_boot, G=G_boot, E=E_boot, pi1=pi1, swap=TRUE, control=con)
   	SPMLE_G_boot[b,] = SPMLE_G_asymp_boot$par
   	SPMLE_E_boot[b,] = SPMLE_E_asymp_boot$par

    ### Fusion (combo) asymptotic (to calculate bootstrap SE even when using asymptotic point estimate)
    ## Define matrices
    Omega_all_boot = c(SPMLE_G_asymp_boot$par, SPMLE_E_asymp_boot$par)
    X_mat = rbind(diag(length_Omega), diag(length_Omega))
    zero_mat = matrix(0, nrow=length_Omega, ncol=length_Omega)

    ## Block diagonal with both inverse hessians
    H_all_boot = -rbind(cbind(SPMLE_G_asymp_boot$H_inv, zero_mat), cbind(zero_mat, SPMLE_E_asymp_boot$H_inv))

    ## Sigma block matrix
    Sigma_GG_boot  = SPMLE_G_asymp_boot$Sigma
    Sigma_EE_boot  = SPMLE_E_asymp_boot$Sigma
    Sigma_GE_boot  = ((ncontrol-1)*cov(SPMLE_G_asymp_boot$zeta0, SPMLE_E_asymp_boot$zeta0) + (ncase-1)*cov(SPMLE_G_asymp_boot$zeta1, SPMLE_E_asymp_boot$zeta1))/n
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
  Omega_all = c(SPMLE_G_par, SPMLE_E_par)
  X_mat = rbind(diag(length_Omega), diag(length_Omega))

  ## Return estimates and Lambda matrix
  return(list(boot_SE = symple_fusion_asymp_boot_SE,
  						SPMLE_G_boot = SPMLE_G_boot,
  						SPMLE_E_boot = SPMLE_E_boot,
  						symple_fusion_boot = symple_fusion_boot
  						))
}


################################################################################ Composite SPMLE

## Sums the neg loglik from SPMLE and swapped likelihoods
composite_lik = function(Omega, D, G, E, pi1) {
  nG = NCOL(G)
  nE = NCOL(E)
  swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
  reverse_swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))
  SPMLE_G_lik = neglikgrad(Omega=Omega, D=D, G=G, E=E, pi1=pi1)  # neg loglik from SPMLE
  SPMLE_E_lik = neglikgrad(Omega=Omega[reverse_swap_order], D=D, G=E, E=G, pi1=pi1)  # neg loglik from swap
  assign(x="grad_val", value=SPMLE_G_lik$gradient + SPMLE_E_lik$gradient[swap_order], pos=parent.frame())  # Store the gradient in the calling environment to retrieve it when ucminf tries to calculate the gradient
  # grad_val <<- SPMLE_G_lik$gradient + SPMLE_E_lik$gradient[swap_order] # Store the gradient in .GlobalEnv to retrieve it when ucminf tries to calculate the gradient
  return(SPMLE_G_lik$objective + SPMLE_E_lik$objective)
}


## Calculate the composite SPMLE (maximize sum of normal SPMLE and swapped likelihoods)
composite_est = function(Omega_start, D, G, E, pi1, control=list()) {
  ## Set control parameters
  con = list(nboot=0, trace=0, usehess=FALSE, max_grad_tol=0.001, num_tries=21)
  con[(names(control))] = control

  nG = NCOL(G)
  nE = NCOL(E)
  swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
  reverse_swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))

  if(con$usehess==TRUE) {
    ## Calculate hessian
    hg = neghess(Omega=Omega_start, D=D, G=G, E=E, pi1=pi1)
    he = neghess(Omega=Omega_start[reverse_swap_order], D=D, G=E, E=G, pi1=pi1)
    H0 = (hg + he[swap_order,swap_order])

    ## Lower triangle of inverse hessian (use identity matrix if hessian is not PD)
    invhessian.lt = tryCatch(solve(H0), error=function(x) diag(length(Omega_start)))[lower.tri(diag(length(Omega_start)),diag=TRUE)]

    ## Maximize composite likelihood
    composite_SPMLE = ucminf::ucminf(par=Omega_start, fn=composite_lik, gr=grad_fn, control=list(trace=con$trace, invhessian.lt=invhessian.lt), D=D, G=G, E=E, pi1=pi1)
  } else {
    composite_SPMLE = ucminf::ucminf(par=Omega_start, fn=composite_lik, gr=grad_fn, control=list(trace=con$trace), D=D, G=G, E=E, pi1=pi1)
  }

  ## Check for convergance, and rerun optimization if necessary
  if(!is.finite(composite_SPMLE$info[1]) | composite_SPMLE$info[1]>con$max_grad_tol) {
    startvals=cbind(rep(0,length(Omega_start)), # create a matrix of starting values
                    matrix(runif(length(Omega_start)*ceiling((con$num_tries-1)/4), min=-1, max=1), nrow=length(Omega_start)) + Omega_start,
                    matrix(rnorm(length(Omega_start)*ceiling((con$num_tries-1)/4)), nrow=length(Omega_start)),
                    matrix(rt(length(Omega_start)*ceiling((con$num_tries-1)/2), df=2), nrow=length(Omega_start)))
    for(j in 1:ncol(startvals)) { # retry optimization
      updatetxt = paste0("UCMINF retry ",j,"\n")
      if(con$trace>-1) {cat(updatetxt)}
      composite_SPMLE = ucminf::ucminf(par=startvals[,j], fn=composite_lik, gr=grad_fn, control=list(trace=con$trace), D=D, G=G, E=E, pi1=pi1)
      if(is.finite(composite_SPMLE$info[1]) & composite_SPMLE$info[1]<con$max_grad_tol) break  # stop once we have convergence
      if(j == ncol(startvals)) stop("UCMINF failed to converge")
    }  # end of j loop
  }
  return(composite_SPMLE)
}


## Calculate asymptotic SE for NPMLE
composite_asymp = function(D, G, E, pi1, Omega_start=NULL, control=list()){
  ## Set control parameters
  con = list(nboot=0, trace=0, usehess=FALSE)
  con[(names(control))] = control

  ## Sizes of arrays
  n = length(D)
  ncase = sum(D)
  ncontrol = n - ncase
  G = as.matrix(G)
  E = as.matrix(E)
  nG = NCOL(G)
  nE = NCOL(E)

  swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
  reverse_swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))

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
  hess_zeta_E = hesszeta(Omega=composite_par[reverse_swap_order], D=D, G=E, E=G, pi1=pi1)

  ## Combined hessian
  H = (hess_zeta_G$hessian + hess_zeta_E$hessian[swap_order,swap_order])/2
  H_inv = -n*solve(H)             # = -(hessian/n)^(-1) = (Gamma_1 - Gamma_2)^-1

  ## Combined zetas
  Zeta0 = (hess_zeta_G$zeta0 + hess_zeta_E$zeta0[, swap_order])/2
  Zeta1 = (hess_zeta_G$zeta1 + hess_zeta_E$zeta1[, swap_order])/2
  Sigma = ((ncontrol-1)*cov(Zeta0) + (ncase-1)*cov(Zeta1))/n

  ## Asymptotic SE for the composite estimator follows the same pattern as the standard SPMLE
  Lambda = H_inv %*% Sigma %*% t(H_inv)  # covar matrix of sqrt(n) * OmegaHat
  SE_asy = sqrt(diag(Lambda)/n)

  ## Return results
  return(list(par   = composite_par,
              SE    = SE_asy,
              cov   = Lambda,
              H_inv = H_inv,
              Sigma = Sigma,
              zeta0 = Zeta0,
              zeta1 = Zeta1,
              H     = H)
  )
}
