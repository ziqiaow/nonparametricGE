########### Functions to provide methods for S3 objects of class spmle

## Return the loglikelihood of spmle objects
logLik.spmle = function(x) x$logLik

## Print spmle objects
print.spmle = function(x, ...) stats:::print.glm(x, ...)

## Summarize spmle objects (similar to summary.glm)
summary.spmle = function (object, correlation = FALSE, symbolic.cor = FALSE, ...) {
  df.r = object$df.residual
  aliased = is.na(coef(object))
  p = object$rank
  p1 = 1L:p
  coef.p = object$coefficients
  covmat = object$cov
  dimnames(covmat) = list(names(coef.p), names(coef.p))
  var.cf = diag(covmat)
  s.err = sqrt(var.cf)
  tvalue = coef.p/s.err
  dn = c("Estimate", "Std. Error")
  pvalue = 2 * pnorm(-abs(tvalue))
  coef.table = cbind(coef.p, s.err, tvalue, pvalue)
  dimnames(coef.table) = list(names(coef.p), c(dn, "z value", "Pr(>|z|)"))
  df.f = object$rank
  iter = object$ucminf$info[c("neval", "maxgradient")]
  retry = list(retries=object$ucminf$retry, start_vals=object$ucminf$start_vals)
  keep = match(c("call", "terms", "family", "deviance", "aic", "contrasts",
                 "df.residual", "null.deviance", "df.null", "na.action"), names(object), 0L)
  ans = c(object[keep], list(deviance.resid = residuals(object, type = "deviance"),
                             coefficients = coef.table,
                             aliased = aliased,
                             df = c(object$rank, df.r, df.f),
                             cov.unscaled = covmat,
                             cov.scaled = covmat,
                             iter = iter,
                             retry = retry))
  if (correlation) {
    dd = sqrt(diag(covmat))
    ans$correlation = covmat/outer(dd, dd)
    ans$symbolic.cor = symbolic.cor
  }
  class(ans) = "summary.spmle"
  return(ans)
}

## Print summary of spmle objects (similar to print.summary.glm)
print.summary.spmle = function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
                                signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Deviance Residuals: \n")
  if (x$df.residual > 5) {
    x$deviance.resid <- setNames(quantile(x$deviance.resid, na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
  }
  xx <- zapsmall(x$deviance.resid, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    df <- if ("df" %in% names(x)) {
      x[["df"]]
    } else {NULL}
    if (!is.null(df) && (nsingular <- df[3L] - df[1L])) {
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
    } else {cat("\nCoefficients:\n")}
    coefs <- x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn <- names(aliased)
      coefs <- matrix(NA, length(aliased), 4L, dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] <- x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }
  cat("\n", apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)), " on", format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
  if (nzchar(mess <- naprint(x$na.action))) {cat("  (", mess, ")\n", sep = "")}
  cat("AIC:", format(x$aic, digits = max(4L, digits + 1L)), "\n")
  cat("UCMINF retries: ", x$retry$retries, ", iterations: ", x$iter[1], ", max gradient at convergence: ", format(x$iter[2], digits=max(5L, digits + 1L)), sep="")
  correl <- x$correlation
  if (!is.null(correl)) {
    p <- NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      } else {
        correl <- format(round(correl, 2L), nsmall = 2L, digits = digits)
        correl[!lower.tri(correl)] <- ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}





