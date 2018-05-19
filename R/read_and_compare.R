#' Function to read saved matrices
#'
#' @param filename file name prefix
#' @param simpath path to find simulations
#' @param matrix matrix to read
#'
#' @return a list of matrices
#' @export
read_matrices = function(filename, simpath="./sims/", matrix="asymp_ests$Sigma_GE") {
	## Create a list of all files
	filePattern = paste(filename, "asy")
	fileList = list.files(path=simpath, pattern=filePattern)

	## Gather the specified matrix from each simulation
	matrixList = lapply(fileList, FUN=function(x) {load(paste0(simpath, x)); eval(parse(text=matrix))})

	return(matrixList)
}



#' Function to read saved simulations
#'
#' @param filename file name prefix
#' @param simpath path to find simulations
#'
#' @return A data frame of simulation results
#' @export
read_sims = function(filename, simpath="./sims/") {
  ## Create a list of all files
  filePattern = paste(filename, "sim")
  fileList = list.files(path=simpath, pattern=filePattern)

  ## Read row & column names
  file_row_names = read.csv(paste0(simpath,fileList[[1]]), header=TRUE, stringsAsFactors=FALSE)[,1]
  file_good_rows = setdiff(1:length(file_row_names), grep("_cov", file_row_names))
  file_row_names = file_row_names[file_good_rows] # drop the covariance matrices
  file_col_names = as.character(read.csv(paste0(simpath,fileList[[1]]), header=FALSE, stringsAsFactors=FALSE)[1,-1])
  name_mat = expand.grid(file_col_names, file_row_names)
  all_names = paste(name_mat[,2], name_mat[,1], sep="_")

  ## Create a data frame of all calls by first applying read.csv, then rbind
  sim_df = do.call("rbind", lapply(fileList, function(x) {c(t(read.csv(paste0(simpath,x), header=TRUE)[file_good_rows,-1]))}))
  colnames(sim_df) = all_names
  return(sim_df)
}


#' Function to compare the performance of the various estimators
#'
#' @param sim_df data frame of simulation results
#' @param truevec true parameter values, with kappa as the intercept.  User must provide truevec OR separate betas.  Default \code{NULL}.
#' @param all_ests Logical: include all estimates. If \code{FALSE}, leave out simple combo, weighted combo, and combo with bootstrapped Lambda. Default \code{TRUE}.
#' @param pi1 population disease rate.  Default \code{NULL}.
#' @param ncase,ncontrol number of cases and controls.  Default \code{NULL}.
#' @param beta0 logistic intercept.  Default \code{NULL}.
#' @param kappa intercept estimated by logistic regression.  If not provided, will be calculated from pi1, ncase, ncontrol, and beta0.  Default \code{NULL}.
#' @param betaG_SNP,betaG_normPRS,betaG_gammaPRS,betaG_bimodalPRS coefficients for genetic variables.  Default \code{NULL}.
#' @param betaE_bin,betaE_norm coefficients for environmental variables.  Default \code{NULL}.
#' @param betaGE_SNP_bin,betaGE_normPRS_bin,betaGE_gammaPRS_bin,betaGE_bimodalPRS_bin,betaGE_SNP_norm,betaGE_normPRS_norm,betaGE_gammaPRS_norm,betaGE_bimodalPRS_norm coefficients for G-E interactions.  Default \code{NULL}.
#'
#' @return a data frame with results of the comparison.  If \code{all_ests=TRUE}, attributes will include a list with coverage matrices for the combination (fusion) estimators.
#' @export
compare_fun = function(sim_df, truevec=NULL, all_ests=TRUE,
											 pi1=NULL, ncase=NULL, ncontrol=NULL, beta0=NULL, kappa=NULL,
											 betaG_SNP=NULL, betaG_normPRS=NULL, betaG_gammaPRS=NULL, betaG_bimodalPRS=NULL,
											 betaE_bin=NULL, betaE_norm=NULL,
											 betaGE_SNP_bin=NULL, betaGE_normPRS_bin=NULL, betaGE_gammaPRS_bin=NULL, betaGE_bimodalPRS_bin=NULL,
											 betaGE_SNP_norm=NULL, betaGE_normPRS_norm=NULL, betaGE_gammaPRS_norm=NULL, betaGE_bimodalPRS_norm=NULL) {

	## Concatenate betas into a single vector
	betas = c(betaG_SNP, betaG_normPRS, betaG_gammaPRS, betaG_bimodalPRS, betaE_bin, betaE_norm,
						betaGE_SNP_bin, betaGE_normPRS_bin, betaGE_gammaPRS_bin, betaGE_bimodalPRS_bin, betaGE_SNP_norm, betaGE_normPRS_norm, betaGE_gammaPRS_norm, betaGE_bimodalPRS_norm)

	## Check that only one of (truevec, betas) were provided
	if(!xor(!is.null(truevec), !is.null(c(pi1, ncase, ncontrol, beta0, kappa, betas)))) stop("Specify either truevec OR betas")

	## Calculate truevec if necessary
	if(is.null(truevec)){
		if(is.null(kappa)) kappa = beta0 + log(ncase*(1-pi1)/(ncontrol*pi1))	# beta0 transformed to kappa because all methods estimate kappa, not beta0
		if(!length(kappa) || !is.finite(kappa)) kappa = 0
		truevec = c(kappa, betas)
	}

  len_v = length(truevec) # length of estimated vector

  # Names of different estimators being compared
  if(any(grepl("boot", colnames(sim_df)))) {  # if the simulation includes bootstraps
    if(all_ests==TRUE) {  # if we want all estimates
      ests_asy_boot_SE = c("logistic", "SPMLE_G", "SPMLE_E", "composite", "simple", "combo_asy")  # methods that have asymptotic & bootstrap estimates of SE
      ests_single_SE = c("weighted", "combo_boot")  # methods that have one estimate of SE
    } else {  # if we don't want simple, weighted, and combo_boot
      ests_asy_boot_SE = c("logistic", "SPMLE_G", "SPMLE_E", "composite", "combo_asy")
      ests_single_SE = NULL
    }
  } else {  # for simulations without bootstraps
    ests_asy_boot_SE = NULL
    if(all_ests==TRUE) {  # if we want all estimates
      ests_single_SE = c("logistic", "SPMLE_G", "SPMLE_E", "composite", "simple", "weighted", "combo_asy")
    } else {  # if we don't want simple and weighted
      ests_single_SE = c("logistic", "SPMLE_G", "SPMLE_E", "composite", "combo_asy")
    }
  }

  # Drop names of estimators that were not used in the simulation being examined
  ests_asy_boot_SE = ests_asy_boot_SE[charmatch(ests_asy_boot_SE, colnames(sim_df), nomatch=-1) == 0]
  ests_single_SE = ests_single_SE[charmatch(ests_single_SE, colnames(sim_df), nomatch=-1) == 0]

  # Names of coefficients estimated
  first_method = c(ests_asy_boot_SE, ests_single_SE)[1]
  coeff_names = substr(colnames(sim_df)[grep(paste0(first_method,"_par_"), colnames(sim_df))], start=nchar(first_method)+6, stop=100)
  names(truevec) = coeff_names

  # Statistics used to compare the estimators
  stats_asy_boot_SE = c("Bias", "SD_of_estimates", "Average_SE_asy", "Average_SE_boot", "MSE", "MSE_efficiency", "CI_length_asy", "CI_length_boot", "Coverage_asy", "Coverage_boot")
  stats_single_SE = c("Bias", "SD_of_estimates", "Average_SE", "MSE", "MSE_efficiency", "CI_length", "Coverage", "cor_SPMLE_G", "cor_SPMLE_E")

  # Create a data frame to hold the comparisons
  comparison = as.data.frame(t(truevec))
  rownames(comparison) = "True_Value"
  coverages = list()

  # MSE of logistic regression - used to measure the performance of other methods
  logistic_MSE=apply(sweep(sim_df[,grep("logistic_par", colnames(sim_df))], 2, STATS=truevec, FUN="-")^2, 2, mean)

  # Compute statistics for the estimators with asymptotic & bootstrap estimates of SE
  for(i in seq_along(ests_asy_boot_SE)){
    # load each method separately
    ests = sim_df[,grep(paste0(ests_asy_boot_SE[i], "_par"), colnames(sim_df))]
    asy_SE = sim_df[,grep(paste0(ests_asy_boot_SE[i], "_asy_SE"), colnames(sim_df))]
    boot_SE = sim_df[,grep(paste0(ests_asy_boot_SE[i], "_boot_SE"), colnames(sim_df))]

    # Calculate Bias, SD of estimates, Average SE, MSE, and MSE efficiency
    Bias = apply(sweep(ests, 2, STATS=truevec, FUN="-"), 2, mean)
    SD_of_estimates = apply(ests, 2, sd)
    Average_SE_asy = colMeans(asy_SE)
    Average_SE_boot = colMeans(boot_SE)
    MSE = apply(sweep(ests, 2, STATS=truevec, FUN="-")^2, 2, mean)
    MSE_efficiency = logistic_MSE/MSE

    # Make 95% confidence intervals
    CI.U.asy = ests+qnorm(0.975)*asy_SE
    CI.L.asy = ests+qnorm(0.025)*asy_SE
    CI_length_asy = colMeans(CI.U.asy-CI.L.asy)

    CI.U.boot = ests+qnorm(0.975)*boot_SE
    CI.L.boot = ests+qnorm(0.025)*boot_SE
    CI_length_boot = colMeans(CI.U.boot-CI.L.boot)

    # Tally the number of CIs that contain the true value
    cover.mat.asy = matrix(ncol=len_v, nrow=nrow(ests))
    cover.mat.boot = matrix(ncol=len_v, nrow=nrow(ests))
    for(j in 1:nrow(ests)) {
      cover.mat.asy[j,] = (truevec<CI.U.asy[j,]) & (truevec>CI.L.asy[j,])
      cover.mat.boot[j,] = (truevec<CI.U.boot[j,]) & (truevec>CI.L.boot[j,])
    }
    if(ests_asy_boot_SE[i] %in% c("combo_asy")) {
    	coverages[[paste0(ests_asy_boot_SE[i], "_asy")]] = cover.mat.asy
    	coverages[[paste0(ests_asy_boot_SE[i], "_boot")]] = cover.mat.boot
    }

    Coverage_asy = 100 * colSums(cover.mat.asy)/nrow(ests)
    Coverage_boot = 100 * colSums(cover.mat.boot)/nrow(ests)

    # Store these statistics in a list
    statlist = lapply(stats_asy_boot_SE, FUN=function(x) get(x, pos=-1))
    names(statlist) = paste(ests_asy_boot_SE[i], stats_asy_boot_SE, sep="_")

    # Include correlation between standard and swapped SPMLE estimates
    if(ests_asy_boot_SE[i]=="SPMLE_E") {
      statlist = c(statlist, list(diag(cor(ests, sim_df[,grep(paste0(ests_asy_boot_SE[i-1], "_par"), colnames(sim_df))]))))
      names(statlist)[length(statlist)] = paste(ests_asy_boot_SE[i], "cor", ests_asy_boot_SE[i-1],sep="_")
    }

    # Remove MSE Efficiency for Logistic Regression
    if(ests_asy_boot_SE[i] == "logistic") statlist["logistic_MSE_efficiency"] = NULL

    statlist = do.call("rbind", statlist)
    colnames(statlist) = colnames(comparison)

    # Add these statistics to the comparison
    comparison = rbind(comparison, statlist)
  }

  # Compute statistics for estimators that have a single estimate of SE
  for(i in seq_along(ests_single_SE)){
    # load each method separately
    ests = sim_df[,grep(paste0(ests_single_SE[i], "_par"), colnames(sim_df))]
    SE = sim_df[,grep(paste0(ests_single_SE[i], "(_asy)?_SE"), colnames(sim_df))]

    # Calculate Bias, SD of estimates, Average SE, MSE, and MSE efficiency
    Bias = apply(sweep(ests, 2, STATS=truevec, FUN="-"), 2, mean)
    SD_of_estimates = apply(ests, 2, sd)
    Average_SE = colMeans(SE)
    MSE = apply(sweep(ests, 2, STATS=truevec, FUN="-")^2, 2, mean)
    MSE_efficiency = logistic_MSE/MSE

    # Make 95% confidence intervals
    CI.U = ests+qnorm(0.975)*SE
    CI.L = ests+qnorm(0.025)*SE
    CI_length = colMeans(CI.U-CI.L)

    # Tally the number of CIs that contain the true value
    cover.mat = matrix(ncol=len_v, nrow=nrow(ests))
    for(j in 1:nrow(ests)) {cover.mat[j,] = (truevec<CI.U[j,]) & (truevec>CI.L[j,])}
    if(ests_single_SE[i] %in% c("combo_asy", "combo_boot")) coverages[[ests_single_SE[i]]] = cover.mat
    Coverage = 100 * colSums(cover.mat)/nrow(ests)

    # Calculate correlation of estimates with standard and swapped SPMLE estimates
    if(ests_single_SE[i]=="SPMLE_G") {
      cor_SPMLE_G = NULL
    } else {
      cor_SPMLE_G = diag(cor(ests, sim_df[,grep("SPMLE_G_par", colnames(sim_df))]))
    }
    if(ests_single_SE[i]=="SPMLE_E") {
      cor_SPMLE_E = NULL
    } else {
      cor_SPMLE_E = diag(cor(ests, sim_df[,grep("SPMLE_E_par", colnames(sim_df))]))
    }

    # Store these statistics in a list
    statlist = lapply(stats_single_SE, FUN=function(x) get(x, pos=-1))
    names(statlist) = paste(ests_single_SE[i], stats_single_SE, sep="_")

    # Remove MSE Efficiency for Logistic Regression
    if(ests_single_SE[i] == "logistic") statlist["logistic_MSE_efficiency"] = NULL


    # Include correlation between asymptotic and bootstrapped Combo estimates
    if(ests_single_SE[i]=="combo_boot") {
      statlist = c(statlist, list(diag(cor(ests, sim_df[,grep(paste0(ests_single_SE[i-1], "_par"), colnames(sim_df))]))))
      names(statlist)[length(statlist)] = paste(ests_single_SE[i], "corr", ests_single_SE[i-1],sep="_")
    }

    statlist = do.call("rbind", statlist)
    colnames(statlist) = colnames(comparison)

    # Add these statistics to the comparison
    comparison = rbind(comparison, statlist)
  }
  
  if(all_ests==TRUE) attr(comparison,"coverages") = coverages  # add coverages for combo (fusion) methods if all_ests=TRUE

  return(comparison)
}
