#' Simulate data and run the analysis
#'
#' @param X simulation number
#' @param filename file name for log and simulations
#' @param simpath,logpath paths for saving files
#' @param control a list of control parameters.
#'     \code{nboot} number of bootstraps.  Default \code{0}.
#'     \code{trace} is a scalar.  If >-1, tracing information is produced.  Default \code{0}.
#'     \code{usehess} logical: precondition optimization with the hessian.  Default \code{TRUE}.
#'     \code{useseed} logical: use the simulation number as the seed.  Default \code{TRUE}.
#' @param pi1 pipulation disease rate
#' @param ncase,ncontrol number of cases and controls
#' @param beta0 logistic intercept
#' @param betaG_SNP,betaG_normPRS,betaG_gammaPRS,betaG_bimodalPRS coefficients for genetic variables
#' @param betaE_bin,betaE_norm coefficients for environmental variables
#' @param betaGE_SNP_bin,betaGE_normPRS_bin,betaGE_gammaPRS_bin,betaGE_bimodalPRS_bin,betaGE_SNP_norm,betaGE_normPRS_norm,betaGE_gammaPRS_norm,betaGE_bimodalPRS_norm coefficients for G-E interactions
#' @param MAF Mean Allele Frequency
#' @param SNP_cor correlation between successive SNPs
#' @param G_normPRS_cor correlation between multivariate normal genetic PRS variables
#' @param E_bin_freq frequency of E=1 for binary environmental variables
#' @param E_norm_cor correlation between multivariate normal environmental variables
#' @param regress_E_bin_on_G_SNP,regress_E_bin_on_G_normPRS,regress_E_bin_on_G_gammaPRS,regress_E_bin_on_G_bimodalPRS,regress_E_norm_on_G_SNP,regress_E_norm_on_G_normPRS,regress_E_norm_on_G_gammaPRS,regress_E_norm_on_G_bimodalPRS  to violate the G-E independence assumption
#'
#' @return A named row vector of estimates and SEs
#' @export
simfun = function(X, filename, simpath, logpath, control=list(),
                  pi1, ncase, ncontrol, beta0,
                  betaG_SNP, betaG_normPRS, betaG_gammaPRS, betaG_bimodalPRS, betaE_bin, betaE_norm,
                  betaGE_SNP_bin, betaGE_normPRS_bin, betaGE_gammaPRS_bin, betaGE_bimodalPRS_bin,
                  betaGE_SNP_norm, betaGE_normPRS_norm, betaGE_gammaPRS_norm, betaGE_bimodalPRS_norm,
                  MAF, SNP_cor, G_normPRS_cor, E_bin_freq, E_norm_cor,
                  regress_E_bin_on_G_SNP, regress_E_bin_on_G_normPRS, regress_E_bin_on_G_gammaPRS, regress_E_bin_on_G_bimodalPRS,
                  regress_E_norm_on_G_SNP, regress_E_norm_on_G_normPRS, regress_E_norm_on_G_gammaPRS, regress_E_norm_on_G_bimodalPRS
) {
  # Check that ncase & ncontrol both > 0
  stopifnot(all(ncase>0, ncontrol>0))

  # Set control parameters
  con = list(nboot=0, trace=0, usehess=TRUE, useseed=TRUE)
  stopifnot(names(control) %in% names(con))
  con[(names(control))] = control

  # Set seed (for parallel sims using mclapply or dopar), and save the RNG state (seed) if it was not set here (used to reproduce doRNG simulations)
  if(con$useseed==TRUE) {set.seed(X)}
  RanSeed = .Random.seed

  # Correct too-small bootstrap sample
  if(con$nboot!=0 && con$nboot<2) {
    warning("To use the bootstrap, nboot must be >= 2.  Skipping bootstrap.")
    con$nboot = 0
  }

  ## Simulate the data (comples: continuous, multimarker, or both)
  dat = simulate_complex(ncase=ncase, ncontrol=ncontrol, beta0=beta0,
                         betaG_SNP=betaG_SNP, betaG_normPRS=betaG_normPRS, betaG_gammaPRS=betaG_gammaPRS, betaG_bimodalPRS=betaG_bimodalPRS,
                         betaE_bin=betaE_bin, betaE_norm=betaE_norm,
                         betaGE_SNP_bin=betaGE_SNP_bin, betaGE_normPRS_bin=betaGE_normPRS_bin, betaGE_gammaPRS_bin=betaGE_gammaPRS_bin, betaGE_bimodalPRS_bin=betaGE_bimodalPRS_bin,
                         betaGE_SNP_norm=betaGE_SNP_norm, betaGE_normPRS_norm=betaGE_normPRS_norm, betaGE_gammaPRS_norm=betaGE_gammaPRS_norm, betaGE_bimodalPRS_norm=betaGE_bimodalPRS_norm,
                         MAF=MAF, SNP_cor=SNP_cor, G_normPRS_cor=G_normPRS_cor, E_bin_freq=E_bin_freq, E_norm_cor=E_norm_cor,
                         regress_E_bin_on_G_SNP=regress_E_bin_on_G_SNP, regress_E_bin_on_G_normPRS=regress_E_bin_on_G_normPRS, regress_E_bin_on_G_gammaPRS=regress_E_bin_on_G_gammaPRS, regress_E_bin_on_G_bimodalPRS=regress_E_bin_on_G_bimodalPRS,
                         regress_E_norm_on_G_SNP=regress_E_norm_on_G_SNP, regress_E_norm_on_G_normPRS=regress_E_norm_on_G_normPRS, regress_E_norm_on_G_gammaPRS=regress_E_norm_on_G_gammaPRS, regress_E_norm_on_G_bimodalPRS=regress_E_norm_on_G_bimodalPRS,
                         control=con
  )

  ## Asymptotic estimates
  asymp_ests = combo_asymp(D=dat$D, G=dat$G, E=dat$E, pi1=pi1, control=con, current_sim=X, logpath=logpath, filename=filename)

  ## Bootstrap estimates
  if(con$nboot>0) boot_ests = combo_boot(D=dat$D, G=dat$G, E=dat$E, pi1=pi1, SPMLE_G_par=asymp_ests$ests['SPMLE_G_par',], SPMLE_E_par=asymp_ests$ests['SPMLE_E_par',], control=con, current_sim=X, logpath=logpath, filename=filename)

  ## All estimates in a single data frame
  if(con$nboot>0) {
    all_estimates = rbind(logistic_par=asymp_ests$ests['logistic_par',],
                          logistic_asy_SE=asymp_ests$ests['logistic_SE',],
                          logistic_boot_SE=boot_ests$ests['logistic_boot_SE',],
                          SPMLE_G_par=asymp_ests$ests['SPMLE_G_par',],
                          SPMLE_G_asy_SE=asymp_ests$ests['SPMLE_G_SE',],
                          SPMLE_G_boot_SE=boot_ests$ests['SPMLE_G_boot_SE',],
                          SPMLE_E_par=asymp_ests$ests['SPMLE_E_par',],
                          SPMLE_E_asy_SE=asymp_ests$ests['SPMLE_E_SE',],
                          SPMLE_E_boot_SE=boot_ests$ests['SPMLE_E_boot_SE',],
                          composite_par=asymp_ests$ests['composite_par',],
                          composite_asy_SE=asymp_ests$ests['composite_SE',],
                          composite_boot_SE=boot_ests$ests['composite_boot_SE',],
                          simple_par=asymp_ests$ests['simple_par',],
                          simple_asy_SE=asymp_ests$ests['simple_SE',],
                          simple_boot_SE=boot_ests$ests['simple_boot_SE',],
                          weighted_par=asymp_ests$ests['weighted_par',],
                          weighted_SE=asymp_ests$ests['weighted_SE',],
                          combo_asy_par=asymp_ests$ests['combo_par',],
                          combo_asy_asy_SE=asymp_ests$ests['combo_SE',],
    											combo_asy_boot_SE=boot_ests$ests['symple_fusion_asymp_boot_SE',],
                          combo_boot_par=boot_ests$ests['combo_boot_par',],
                          combo_boot_SE=boot_ests$ests['combo_boot_SE',]
    )
  } else {
    all_estimates = rbind(logistic_par=asymp_ests$ests['logistic_par',],
                          logistic_asy_SE=asymp_ests$ests['logistic_SE',],
                          SPMLE_G_par=asymp_ests$ests['SPMLE_G_par',],
                          SPMLE_G_asy_SE=asymp_ests$ests['SPMLE_G_SE',],
                          SPMLE_E_par=asymp_ests$ests['SPMLE_E_par',],
                          SPMLE_E_asy_SE=asymp_ests$ests['SPMLE_E_SE',],
                          composite_par=asymp_ests$ests['composite_par',],
                          composite_asy_SE=asymp_ests$ests['composite_SE',],
                          simple_par=asymp_ests$ests['simple_par',],
                          simple_asy_SE=asymp_ests$ests['simple_SE',],
                          weighted_par=asymp_ests$ests['weighted_par',],
                          weighted_SE=asymp_ests$ests['weighted_SE',],
                          combo_asy_par=asymp_ests$ests['combo_par',],
                          combo_asy_asy_SE=asymp_ests$ests['combo_SE',]
    )
  }

  ## Set variable names for returning as a single row
  onerow = c(t(all_estimates))
  name_mat=expand.grid(c("Intercept", colnames(all_estimates)[-1]), rownames(all_estimates))
  all_names=paste(name_mat[,2], name_mat[,1], sep="_")
  names(onerow) = all_names

  # Write final estimates for this simulation to a file
  write.csv(all_estimates, file=paste0(simpath, filename, " sim", X, ".csv"))

  # Write estimates of Lambda_all for this simulation to a file
  if(con$nboot>0) {
    # write.table(rbind(boot_ests$Lambda_all, asymp_ests$Lambda_all, cbind(asymp_ests$composite_H_inv, asymp_ests$composite_Sigma)),
    #             file=paste0(simpath, filename, " Lambda", X, ".csv"),
    #             row.names=c(rep("boot", 2*ncol(all_estimates)),rep("asy", 2*ncol(all_estimates)),rep("composite_H_inv_&_Sigma", ncol(all_estimates))),
    #             col.names=FALSE, qmethod='double', sep=',')
  	save(asymp_ests, boot_ests, RanSeed, file=paste0(simpath, filename, " asyboot", X, ".Rdata"))
  } else {
    # write.table(rbind(asymp_ests$Lambda_all, cbind(asymp_ests$composite_H_inv, asymp_ests$composite_Sigma)),
    #             file=paste0(simpath, filename, " Lambda", X, ".csv"),
    #             row.names=c(rep("asy", 2*ncol(all_estimates)),rep("composite_H_inv_&_Sigma", ncol(all_estimates))),
    #             col.names=FALSE, qmethod='double', sep=',')
  	save(asymp_ests, RanSeed, file=paste0(simpath, filename, " asy", X, ".Rdata"))
  }
  # Return the estimates as a single vector
  return(onerow)
}
