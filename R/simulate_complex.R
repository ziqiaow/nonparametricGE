#' Simulate case-control data with multivariate, possibly dependent genetic and environmental components.
#'
#' \code{simulate_complex} simulates data to be analyzed by \code{symple}, logistic regression, or other models.
#'
#' This function was used to simulate data for Stalder et. al. (2017).
#' The user can specify up to four types of genetic variables (SNPs with additive effects under HWE,
#' polygenic risk scores with normal, gamma, and bimodal distributions), each of which can be multivariate.
#' Two types of environmental variables (binary and normal) can also be potentially multivariate.
#'
#' If both G and E are multivariate, \code{beta_GE_} and \code{regress_E_} arguments iterate G quickly and E slowly.
#' For example, if G_SNP and E_bin are both bivariate, \code{betaGE_SNP_bin} is ordered (G1*E1, G2*E1, G1*E2, G2*E2).
#'
#' @param ncase,ncontrol number of cases and controls.
#' @param beta0 logistic intercept. Can be manipulated to change the disease rate.
#' @param betaG_SNP,betaG_normPRS,betaG_gammaPRS,betaG_bimodalPRS coefficients for genetic variable main effects.
#'  Genetic variables can include SNPs and polygenic risk scores with normal, gamma, and bimodal distributions.
#'  A value of \code{NULL} for coefficient of a genetic variable means that type of genetic variable will not be generated.
#'  Vector valued coefficients produce multivariate genetic data.
#' @param betaE_bin,betaE_norm coefficients for environmental variable main effects. Environmental variables can include
#'  binary and normally distributed random variables. A value of \code{NULL} for coefficient of an environmental variable
#'  means that type of environmental variable will not be generated. Vector valued coefficients produce multivariate environmental data.
#' @param betaGE_SNP_bin,betaGE_normPRS_bin,betaGE_gammaPRS_bin,betaGE_bimodalPRS_bin,betaGE_SNP_norm,betaGE_normPRS_norm,betaGE_gammaPRS_norm,betaGE_bimodalPRS_norm coefficients for
#'  multiplicative G*E effects. The length of the coefficient of any given G*E interaction must be equal to the product of the lengths
#'  of the coefficients of the corresponding G and E main effects. If one of the main effects is \code{NULL} the interaction must be \code{NULL} as well.
#' @param MAF Minor Allele Frequency of SNPs. This vector is the same length as \code{beta_G_SNP} and has values between 0 and 1.
#'  The MAF is used to generate SNP data that is in Hardy-Weinberg Equilibrium.
#' @param SNP_cor scalar between values between -1 and 1. SNPs are simulated by generating multivariate normal random draws with an AR(\code{SNP_cor})
#'  covariance matrix. These normal draws are then trichotomized according to HWE to simulate SNPs.
#' @param G_normPRS_cor correlation matrix for multivariate normal polygenic risk scores. In the bivariate case, a 2x2 matrix or a scalar
#'  (for correlation) are accepted. NULL values indicate independent polygenic risk scores.
#' @param E_bin_freq marginal probability that E_bin = 1. Must have length equal to \code{length(betaE_bin)} and values between 0 and 1.
#' @param E_norm_cor correlation matrix for multivariate normal environmental variable. In the bivariate case, a 2x2 matrix or a scalar
#'  (for correlation) are accepted. NULL values indicate independent normal environmental variables.
#' @param regress_E_bin_on_G_SNP,regress_E_bin_on_G_normPRS,regress_E_bin_on_G_gammaPRS,regress_E_bin_on_G_bimodalPRS allow the simulation of case-control data that violates
#'  the G-E independence assumption. If these arguments are \code{NULL} or all \code{0}s, the binary environmental variables will be independent
#'  of the genetic variables. If non-null, the length of the regression argument must be equal to the product of the lengths of the coefficients
#'  of the corresponding G and E main effects. The conditional expectations of binary environmental variables depend on the product of the
#'  \code{regress_E_bin} arguments and the corresponding genetic variables.
#' @param regress_E_norm_on_G_SNP,regress_E_norm_on_G_normPRS,regress_E_norm_on_G_gammaPRS,regress_E_norm_on_G_bimodalPRS allow the simulation of case-control data that violates
#'  the G-E independence assumption. If these arguments are \code{NULL} or all \code{0}s, the normal environmental variables will be independent
#'  of the genetic variables. If non-null, the length of the regression argument must be equal to the product of the lengths of the coefficients
#'  of the corresponding G and E main effects. The conditional means of normal environmental variables depend on the product of the
#'  \code{regress_E_norm} arguments and the corresponding genetic variables.
#' @param control a list of control parameters.
#'     \code{trace} is a scalar.  If >-1, tracing information is produced.  Default \code{trace=0}.
#' @return \code{simulate_complex} produces a list with three elements:
#'  \item{\code{D}}{a binary vector with \code{ncontrol} \code{0}s and \code{ncase} \code{1}s.}
#'  \item{\code{G}}{a matrix with \code{ncontrol + ncase} rows and a column for each genetic variable. Genetic variables are ordered: SNPs, normal PRSs, gamma PRSs.}
#'  \item{\code{E}}{a matrix with \code{ncontrol + ncase} rows and a column for each environmental variable. Environmental variables are ordered: binary, normal.}
#' @seealso \code{symple}
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
#' # Simulation with 5 SNPs and a single normal environmental variable
#' # that is dependent on G1 with an R^2 of 0.001.
#' # True population disease rate in this simulation is 0.03.
#' # This simulation scenario was used in the Supplementary Material of Stalder et. al. (2017)
#' dat = simulate_complex(ncase=1000,
#'                        ncontrol=1000,
#'                        beta0=-4.19,
#'                        betaG_SNP=c(log(1.2), log(1.2), 0, log(1.2), 0),
#'                        betaE_norm=(qnorm(0.75)-qnorm(0.25))*c(log(1.5)),
#'                        betaGE_SNP_norm=(qnorm(0.75)-qnorm(0.25))*c(log(1.3), 0, 0, log(1.3), 0),
#'                        MAF=c(0.1, 0.3, 0.3, 0.3, 0.1),
#'                        SNP_cor=0.7,
#'                        regress_E_norm_on_G_SNP=c(sqrt(0.001),rep(0,4)),
#'                        control=list(trace=1))
#' @export
simulate_complex = function(ncase, ncontrol, beta0,
                            betaG_SNP=NULL, betaG_normPRS=NULL, betaG_gammaPRS=NULL, betaG_bimodalPRS=NULL,
                            betaE_bin=NULL, betaE_norm=NULL,
                            betaGE_SNP_bin=NULL, betaGE_normPRS_bin=NULL, betaGE_gammaPRS_bin=NULL, betaGE_bimodalPRS_bin=NULL,
                            betaGE_SNP_norm=NULL, betaGE_normPRS_norm=NULL, betaGE_gammaPRS_norm=NULL, betaGE_bimodalPRS_norm=NULL,
                            MAF=NULL, SNP_cor=0, G_normPRS_cor=0, E_bin_freq=0.5, E_norm_cor=0,
                            regress_E_bin_on_G_SNP=NULL, regress_E_bin_on_G_normPRS=NULL, regress_E_bin_on_G_gammaPRS=NULL, regress_E_bin_on_G_bimodalPRS=NULL,
                            regress_E_norm_on_G_SNP=NULL, regress_E_norm_on_G_normPRS=NULL, regress_E_norm_on_G_gammaPRS=NULL, regress_E_norm_on_G_bimodalPRS=NULL,
                            control=list()) {
  # Set control parameters
  con = list(trace=get0("trace", envir=as.environment(control), ifnotfound=0))

  ### Set the number of patients & variables
  ncontrol      = ceiling(ncontrol) # force ncontrol to be an integer
  ncase         = ceiling(ncase) # force ncase to be an integer
  nG_SNP        = length(betaG_SNP)
  nG_normPRS    = length(betaG_normPRS)
  nG_gammaPRS   = length(betaG_gammaPRS)
  nG_bimodalPRS = length(betaG_bimodalPRS)
  nE_bin        = length(betaE_bin)
  nE_norm       = length(betaE_norm)
  G_normPRS_cor = as.matrix(G_normPRS_cor)
  E_norm_cor    = as.matrix(E_norm_cor)
  tol           = sqrt(.Machine$double.eps)  # set tolerance (used for checking for positive definiteness, etc)

  #### Check for errors ####
  if(ncontrol<1 | ncase<1) stop("Must request at least 1 case and 1 control")
  if(sum(nG_SNP,nG_normPRS,nG_gammaPRS,nG_bimodalPRS)==0 | sum(nE_bin,nE_norm)==0)   stop("Must include at least 1 Genetic and 1 Environmental variable")
  if(length(betaGE_SNP_bin)         != length(betaG_SNP)*length(betaE_bin))          stop("length(betaGE_SNP_bin) must equal length(betaG_SNP)*length(betaE_bin)")
  if(length(betaGE_normPRS_bin)     != length(betaG_normPRS)*length(betaE_bin))      stop("length(betaGE_normPRS_bin) must equal length(betaG_normPRS)*length(betaE_bin)")
  if(length(betaGE_gammaPRS_bin)    != length(betaG_gammaPRS)*length(betaE_bin))     stop("length(betaGE_gammaPRS_bin) must equal length(betaG_gammaPRS)*length(betaE_bin)")
  if(length(betaGE_bimodalPRS_bin)  != length(betaG_bimodalPRS)*length(betaE_bin))   stop("length(betaGE_bimodalPRS_bin) must equal length(betaG_bimodalPRS)*length(betaE_bin)")
  if(length(betaGE_SNP_norm)        != length(betaG_SNP)*length(betaE_norm))         stop("length(betaGE_SNP_norm) must equal length(betaG_SNP)*length(betaE_norm)")
  if(length(betaGE_normPRS_norm)    != length(betaG_normPRS)*length(betaE_norm))     stop("length(betaGE_normPRS_norm) must equal length(betaG_normPRS)*length(betaE_norm)")
  if(length(betaGE_gammaPRS_norm)   != length(betaG_gammaPRS)*length(betaE_norm))    stop("length(betaGE_gammaPRS_norm) must equal length(betaG_gammaPRS)*length(betaE_norm)")
  if(length(betaGE_bimodalPRS_norm) != length(betaG_bimodalPRS)*length(betaE_norm))  stop("length(betaGE_bimodalPRS_norm) must equal length(betaG_bimodalPRS)*length(betaE_norm)")
  if(!is.null(regress_E_bin_on_G_SNP))         if(length(regress_E_bin_on_G_SNP)         != length(betaG_SNP)*length(betaE_bin))          stop("length(regress_E_bin_on_G_SNP) must equal length(betaG_SNP)*length(betaE_bin)")
  if(!is.null(regress_E_bin_on_G_normPRS))     if(length(regress_E_bin_on_G_normPRS)     != length(betaG_normPRS)*length(betaE_bin))      stop("length(regress_E_bin_on_G_normPRS) must equal length(betaG_normPRS)*length(betaE_bin)")
  if(!is.null(regress_E_bin_on_G_gammaPRS))    if(length(regress_E_bin_on_G_gammaPRS)    != length(betaG_gammaPRS)*length(betaE_bin))     stop("length(regress_E_bin_on_G_gammaPRS) must equal length(betaG_gammaPRS)*length(betaE_bin)")
  if(!is.null(regress_E_bin_on_G_bimodalPRS))  if(length(regress_E_bin_on_G_bimodalPRS)  != length(betaG_bimodalPRS)*length(betaE_bin))   stop("length(regress_E_bin_on_G_bimodalPRS) must equal length(betaG_bimodalPRS)*length(betaE_bin)")
  if(!is.null(regress_E_norm_on_G_SNP))        if(length(regress_E_norm_on_G_SNP)        != length(betaG_SNP)*length(betaE_norm))         stop("length(regress_E_norm_on_G_SNP) must equal length(betaG_SNP)*length(betaE_norm)")
  if(!is.null(regress_E_norm_on_G_normPRS))    if(length(regress_E_norm_on_G_normPRS)    != length(betaG_normPRS)*length(betaE_norm))     stop("length(regress_E_norm_on_G_normPRS) must equal length(betaG_normPRS)*length(betaE_norm)")
  if(!is.null(regress_E_norm_on_G_gammaPRS))   if(length(regress_E_norm_on_G_gammaPRS)   != length(betaG_gammaPRS)*length(betaE_norm))    stop("length(regress_E_norm_on_G_gammaPRS) must equal length(betaG_gammaPRS)*length(betaE_norm)")
  if(!is.null(regress_E_norm_on_G_bimodalPRS)) if(length(regress_E_norm_on_G_bimodalPRS) != length(betaG_bimodalPRS)*length(betaE_norm))  stop("length(regress_E_norm_on_G_bimodalPRS) must equal length(betaG_bimodalPRS)*length(betaE_norm)")
  if(length(MAF) != length(betaG_SNP)) stop("length(MAF) must equal length(betaG_SNP)")
  if(!is.null(MAF)) if(max(MAF)>=1 | min(MAF)<=0) stop("All values of MAF must be between 0 and 1")
  if(!is.null(betaG_SNP)) if(!all(length(SNP_cor)==1 & SNP_cor<1 & SNP_cor>-1)) stop("SNP_cor must be between -1 and 1")
  if(!is.null(betaE_bin)) {
  	if(!all(E_bin_freq<1 & E_bin_freq>0)) {stop("All E_bin_freq values must be between 0 and 1")}
  	if(length(E_bin_freq) != length(betaE_bin)) {
  		if(length(E_bin_freq)==1) E_bin_freq = rep(E_bin_freq, times=length(betaE_bin))  # special case: multiple binary environmental variables with the same frequency
  		else stop("length(E_bin_freq) must equal length(betaE_bin) or 1")
  	}
  }

  ### Create the beta (coefficient) vectors
  ## First the main effects
  betaG = c(betaG_SNP, betaG_normPRS, betaG_gammaPRS, betaG_bimodalPRS)
  betaE = c(betaE_bin, betaE_norm)
  ## Then the interaction, iterating G faster than E
  betaGE = NULL
  S = Pn = Pg = Pb = 0
  for(i in seq_along(betaE_bin)) {
    for(j in seq_along(betaG_SNP)) {
      S = S+1
      betaGE = c(betaGE, betaGE_SNP_bin[S])
    }
    for(j in seq_along(betaG_normPRS)) {
      Pn = Pn+1
      betaGE = c(betaGE, betaGE_normPRS_bin[Pn])
    }
    for(j in seq_along(betaG_gammaPRS)) {
      Pg = Pg+1
      betaGE = c(betaGE, betaGE_gammaPRS_bin[Pg])
    }
    for(j in seq_along(betaG_bimodalPRS)) {
      Pb = Pb+1
      betaGE = c(betaGE, betaGE_bimodalPRS_bin[Pb])
    }
  }
  S = Pn = Pg = Pb = 0
  for(i in seq_along(betaE_norm)) {
    for(j in seq_along(betaG_SNP)) {
      S = S+1
      betaGE = c(betaGE, betaGE_SNP_norm[S])
    }
    for(j in seq_along(betaG_normPRS)) {
      Pn = Pn+1
      betaGE = c(betaGE, betaGE_normPRS_norm[Pn])
    }
    for(j in seq_along(betaG_gammaPRS)) {
      Pg = Pg+1
      betaGE = c(betaGE, betaGE_gammaPRS_norm[Pg])
    }
    for(j in seq_along(betaG_bimodalPRS)) {
      Pb = Pb+1
      betaGE = c(betaGE, betaGE_bimodalPRS_norm[Pb])
    }
  }

  ##### Create covariance matrix and allele frequency matrix for SNP genetic variables #####
  chol_Sigma_G_SNP = chol(SNP_cor^abs(outer(1:nG_SNP, 1:nG_SNP, "-")))  # Cholesky decomposition of covariance matrix for G variables
  HWE = apply(rbind((1-MAF)^2, 2*MAF*(1-MAF)), 2, cumsum)  # Cumulative probability of genotypes aa and aA under Hardy-Weinberg equilibrium

  ### create the cholesky decomposition of the correlation matrix of the normal PRS variables
  if(length(G_normPRS_cor)==1 & length(betaG_normPRS)==2) {  # special case: 2 normal E variables and a scalar correlation
    chol_Sigma_G_normPRS = chol(matrix(c(1,G_normPRS_cor,G_normPRS_cor,1), ncol=2, nrow=2))
  } else if(dim(G_normPRS_cor)[1]==length(betaG_normPRS)  # is the given correlation matrix the correct dimension
            & ifelse(test=isSymmetric(G_normPRS_cor), yes=all(eigen(G_normPRS_cor)$values+tol > 0), no=FALSE)  # is it positive definite
            & all(abs(diag(G_normPRS_cor)-rep(1, times=dim(G_normPRS_cor)[1])) <= tol)  # are the diagonal elements all 1
            & all(abs(G_normPRS_cor)-tol <= 1) ) {  # are all elements between -1 and 1
    chol_Sigma_G_normPRS = chol(G_normPRS_cor)
  } else {  # if G_normPRS_cor is not a proper correlation matrix, use the identity matrix instead
    chol_Sigma_G_normPRS = diag(length(betaG_normPRS))
    if(!(is.null(betaG_normPRS) | G_normPRS_cor==0)) warning("Incorrect correlation matrix for normal PRS variables: assuming no correlation")
  }

  ### create the cholesky decomposition of the correlation matrix of the normal environmental variables
  if(length(E_norm_cor)==1 & length(betaE_norm)==2) {  # special case: 2 normal E variables and a scalar correlation
    chol_Sigma_E = chol(matrix(c(1,E_norm_cor,E_norm_cor,1), ncol=2, nrow=2))
  } else if(dim(E_norm_cor)[1]==length(betaE_norm)  # is the given correlation matrix the correct dimension
            & ifelse(isSymmetric(E_norm_cor), all(eigen(E_norm_cor)$values+tol > 0), FALSE)  # is it positive definite
            & all(abs(diag(E_norm_cor)-rep(1, times=dim(E_norm_cor)[1])) <= tol)  # are the diagonal elements all 1
            & all(abs(E_norm_cor)-tol <= 1) ) {  # are all elements between -1 and 1
    chol_Sigma_E = chol(E_norm_cor)
  } else {  # if E_norm_cor is not a proper correlation matrix, use the identity matrix instead
    chol_Sigma_E = diag(length(betaE_norm))
    if(!(is.null(betaE_norm) | E_norm_cor==0)) warning("Incorrect correlation matrix for normal environmental variables: assuming no correlation")
  }

  ##### Initialize variables for simulation #####
  ntotal   = ncontrol + ncase  # starting size of cohort, from which sample is drawn
  p_hat    = 0.05  # population disease rate of 5% to start
  D        = 0  # Initialize D

  ##### Generate the cohort & select the case-control sample from it
  while(length(D)<ncontrol+ncase) {
    ntotal  = ceiling(ntotal*(1+(ncase/(p_hat*ntotal))))  # increase cohort size until we get the desired sample size

    ## Simulate genetic variables
    if(nG_SNP>0) {
      G_sim_SNP = matrix(rnorm(ntotal*nG_SNP, mean=0, sd=1), ncol=nG_SNP) %*% chol_Sigma_G_SNP  # simulating SNPs - first as a multivariate normal variable
      G_sim_SNP = apply(rbind(HWE, G_sim_SNP), 2, FUN=function(x) cut(x[-(1:2)], breaks=c(-Inf,qnorm(x[1:2]), Inf), labels=FALSE)-1)  # trichotomize, following HWE
    } else G_sim_SNP = NULL
    if(nG_normPRS>0) {  # simulating PRS G as N(0,1)
      G_sim_normPRS = (matrix(rnorm(ntotal*nG_normPRS, mean=0, sd=1),ncol=nG_normPRS) %*% chol_Sigma_G_normPRS)
    } else {G_sim_normPRS = NULL}
    if(nG_gammaPRS>0) {  # simulating PRS G from a Gamma(shape=20, rate=20) which has mean = 1
      G_sim_gammaPRS = matrix(rgamma(ntotal*nG_gammaPRS,20,20),ncol=nG_gammaPRS)
    } else {G_sim_gammaPRS = NULL}
    if(nG_bimodalPRS>0) {  # simulating PRS G from a bimodal mixture of normals which has mean = 0 and SD = 1
      G_sim_bimodalPRS = matrix(sample(c(rnorm(n=(ntotal*nG_bimodalPRS+1)/2, mean=-0.85, sd=0.55), rnorm(n=(ntotal*nG_bimodalPRS+1)/2, mean=0.85, sd=0.55)), size=ntotal*nG_bimodalPRS), ncol=nG_bimodalPRS)
    } else {G_sim_bimodalPRS = NULL}
    G_sim   = cbind(G_sim_SNP, G_sim_normPRS, G_sim_gammaPRS, G_sim_bimodalPRS)

    ## Simulate environmental variables
    E_sim_bin = E_sim_norm = NULL
    if(nE_bin>0) {  # simulate binary environmental variables, if there are any
      if(is.null(regress_E_bin_on_G_SNP) & is.null(regress_E_bin_on_G_normPRS) & is.null(regress_E_bin_on_G_gammaPRS) & is.null(regress_E_bin_on_G_bimodalPRS)) {  # if binary E is independent of G, proceed normally
        for(i in 1:nE_bin) {
          E_sim_bin = as.matrix(cbind(E_sim_bin,rbinom(n=ntotal, size=1, prob=E_bin_freq[i])))  # simulating binary E
        }
      } else {  # If G and E are not independent, include correlation
        for(i in 1:nE_bin) {
          GE_interact = numeric(ntotal)  # Vector to store the effect of the G*E correlation
          if(!is.null(regress_E_bin_on_G_SNP)) {  # If binary E isn't independent of SNPs, include that correlation
            for(j in 1:nG_SNP) {
              GE_interact = GE_interact + scale(G_sim_SNP[,j]) * regress_E_bin_on_G_SNP[(i-1)*nG_SNP+j]
            }
          }
          if(!is.null(regress_E_bin_on_G_normPRS)) {  # If binary E isn't independent of normal G, include that correlation
            for(j in 1:nG_normPRS) {
              GE_interact = GE_interact + scale(G_sim_normPRS[,j]) * regress_E_bin_on_G_normPRS[(i-1)*nG_normPRS+j]
            }
          }
          if(!is.null(regress_E_bin_on_G_gammaPRS)) {  # If binary E isn't independent of gamma G, include that correlation
            for(j in 1:nG_gammaPRS) {
              GE_interact = GE_interact + scale(G_sim_gammaPRS[,j]) * regress_E_bin_on_G_gammaPRS[(i-1)*nG_gammaPRS+j]
            }
          }
          if(!is.null(regress_E_bin_on_G_bimodalPRS)) {  # If binary E isn't independent of bimodal G, include that correlation
            for(j in 1:nG_bimodalPRS) {
              GE_interact = GE_interact + scale(G_sim_bimodalPRS[,j]) * regress_E_bin_on_G_bimodalPRS[(i-1)*nG_bimodalPRS+j]
            }
          }
          prob_E_bin = plogis(GE_interact + qlogis(E_bin_freq[i]))  # Vector of probabilities for binary E
          E_sim_bin = as.matrix(cbind(E_sim_bin,rbinom(n=ntotal, size=1, prob=prob_E_bin)))  # simulating binary E that is correlated with G
        }
      }
    }

    if(nE_norm>0) {
      E_sim_norm = matrix(rnorm(ntotal*nE_norm, mean=0, sd=1),ncol=nE_norm) %*% chol_Sigma_E  # simulating continuous E from a N(0,1) distribution with correlation E_norm_cor
      # If normal E isn't independent of G, add that correlation
      for(i in 1:nE_norm) {
        # GE_interact = numeric(ntotal)  # Vector to store the effect of the G*E correlation
        if(!is.null(regress_E_norm_on_G_SNP)) {  # If normal E isn't independent of SNPs, include that correlation
          for(j in 1:nG_SNP) {
            E_sim_norm[,i] = E_sim_norm[,i] + scale(G_sim_SNP[,j]) * regress_E_norm_on_G_SNP[(i-1)*nG_SNP+j]
          }
        }
        if(!is.null(regress_E_norm_on_G_normPRS)) {  # If normal E isn't independent of normal G, include that correlation
          for(j in 1:nG_normPRS) {
            E_sim_norm[,i] = E_sim_norm[,i] + scale(G_sim_normPRS[,j]) * regress_E_norm_on_G_normPRS[(i-1)*nG_normPRS+j]
          }
        }
        if(!is.null(regress_E_norm_on_G_gammaPRS)) {  # If normal E isn't independent of gamma G, include that correlation
          for(j in 1:nG_gammaPRS) {
            E_sim_norm[,i] = E_sim_norm[,i] + scale(G_sim_gammaPRS[,j]) * regress_E_norm_on_G_gammaPRS[(i-1)*nG_gammaPRS+j]
          }
        }
        if(!is.null(regress_E_norm_on_G_bimodalPRS)) {  # If normal E isn't independent of bimodal G, include that correlation
          for(j in 1:nG_bimodalPRS) {
            E_sim_norm[,i] = E_sim_norm[,i] + scale(G_sim_bimodalPRS[,j]) * regress_E_norm_on_G_bimodalPRS[(i-1)*nG_bimodalPRS+j]
          }
        }
      }
    }

    E_sim   = cbind(E_sim_bin,E_sim_norm)

    prob    = plogis((model.matrix( ~ G_sim*E_sim) %*% c(beta0,betaG,betaE,betaGE)))  # inverse logit of linear combination of variables & betas
    p_hat   = max(mean(prob), 0.0001)  # average probability of D=1 in cohort. Min of 0.0001 prevents division by 0 when simulating very rare diseases
    D_sim   = rbinom(ntotal,1,prob)  # 1=case, 0=control
    D       = c(head(D_sim[D_sim==0],ncontrol),head(D_sim[D_sim==1],ncase))  # select the first ncontrol controls and ncase cases, if available
  }  # End of while loop

  ## select the first ncase cases and ncontrol controls
  G = rbind(as.matrix(head(G_sim[D_sim==0,],ncontrol)),as.matrix(head(G_sim[D_sim==1,],ncase)))
  E = rbind(as.matrix(head(E_sim[D_sim==0,],ncontrol)),as.matrix(head(E_sim[D_sim==1,],ncase)))

  ## Print details about the simulation if requested
  if(con$trace > -1) {
    # Print the regression of E on G, if they are not independent
    if(!is.null(regress_E_bin_on_G_SNP)) {
      updatetxt=paste0("\nregress_E_bin_on_G_SNP = ", paste(regress_E_bin_on_G_SNP, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(glm(E_sim_bin ~ G_sim_SNP, family=binomial)))
    }
    if(!is.null(regress_E_bin_on_G_normPRS)) {
      updatetxt=paste0("\nregress_E_bin_on_G_normPRS = ", paste(regress_E_bin_on_G_normPRS, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(glm(E_sim_bin ~ G_sim_normPRS, family=binomial)))
    }
    if(!is.null(regress_E_bin_on_G_gammaPRS)) {
      updatetxt=paste0("\nregress_E_bin_on_G_gammaPRS = ", paste(regress_E_bin_on_G_gammaPRS, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(glm(E_sim_bin ~ G_sim_gammaPRS, family=binomial)))
    }
    if(!is.null(regress_E_bin_on_G_bimodalPRS)) {
      updatetxt=paste0("\nregress_E_bin_on_G_bimodalPRS = ", paste(regress_E_bin_on_G_bimodalPRS, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(glm(E_sim_bin ~ G_sim_bimodalPRS, family=binomial)))
    }
    if(!is.null(regress_E_norm_on_G_SNP)) {
      updatetxt=paste0("\nregress_E_norm_on_G_SNP = ", paste(regress_E_norm_on_G_SNP, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(lm(E_sim_norm ~ G_sim_SNP)))
    }
    if(!is.null(regress_E_norm_on_G_normPRS)) {
      updatetxt=paste0("\nregress_E_norm_on_G_normPRS = ", paste(regress_E_norm_on_G_normPRS, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(lm(E_sim_norm ~ G_sim_normPRS)))
    }
    if(!is.null(regress_E_norm_on_G_gammaPRS)) {
      updatetxt=paste0("\nregress_E_norm_on_G_gammaPRS = ", paste(regress_E_norm_on_G_gammaPRS, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(lm(E_sim_norm ~ G_sim_gammaPRS)))
    }
    if(!is.null(regress_E_norm_on_G_bimodalPRS)) {
      updatetxt=paste0("\nregress_E_norm_on_G_bimodalPRS = ", paste(regress_E_norm_on_G_bimodalPRS, collapse=" "), "\n")
      cat(updatetxt)
      print(summary(lm(E_sim_norm ~ G_sim_bimodalPRS)))
    }

    # Print variable correlations when applicable
    if(nG_SNP>1) {
      updatetxt=paste("\nCorrelation of SNPs:", paste(cor(G_sim_SNP)[cbind(1:(nG_SNP-1), 2:nG_SNP)], collapse=" "), "\n")
      cat(updatetxt)
    }
    if(nG_normPRS>1) {
      if(nG_normPRS==2) {
        updatetxt = paste("\nCorrelation of normal PRS:", cor(G_sim_normPRS[,1],G_sim_normPRS[,2]), "\n")
        cat(updatetxt)
      }
      else {
        cat("\nCorrelation of normal PRS:\n")
        print(cor(G_sim_normPRS))
      }
    }
    if(nE_norm>1) {
      if(nE_norm==2) {
        updatetxt = paste("\nCorrelation of normal E variables:", cor(E_sim_norm[,1],E_sim_norm[,2]), "\n")
        cat(updatetxt)
      }
      else {
        cat("\nCorrelation of normal E variables:\n")
        print(cor(E_sim_norm))
      }
    }

    # Print disease rate
    updatetxt = paste("\nDisease prevalance:", mean(D_sim), "\n\n")
    cat(updatetxt)
  }

  #### Bind D, G, and E into a list
  dat = list(D=D, G=G, E=E)
  return(dat)
}
