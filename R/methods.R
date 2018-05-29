########### Functions to provide methods for S3 objects of class spmle


## Return the loglikelihood of spmle objects
#' @export
logLik.spmle = function(object, ...) object$logLik


## Print spmle objects
#' @export
print.spmle = function(x, digits=max(3L, getOption("digits") - 3L), ...) {
  cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  if(length(coef(x))) {
    cat("Coefficients")
    if(is.character(co <- x$contrasts)) {
      cat("  [contrasts: ", apply(cbind(names(co), co), 1L, paste, collapse = "="), "]")
    }
    cat(":\n")
    print.default(format(x$coefficients, digits = digits), print.gap = 2, quote = FALSE)
  } else {cat("No coefficients\n\n")}
  cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ", x$df.residual, "Residual\n")
  if(nzchar(mess <- naprint(x$na.action))) {cat("  (", mess, ")\n", sep = "")}
  if(!is.null(x$deviance)) {
    cat("Null Deviance:\t   ", format(signif(x$null.deviance, digits)), "\nResidual Deviance:", format(signif(x$deviance, digits)), "\tAIC:", format(signif(x$aic, digits)))
    cat("\n")
  }
  invisible(x)
}


## Summarize spmle objects (similar to summary.glm)
#' @export
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
  ans = c(object[keep], list(pearson.resid = residuals(object),
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
#' @export
print.summary.spmle = function (x, digits = max(3L, getOption("digits") - 3L), symbolic.cor = x$symbolic.cor,
                                signif.stars = getOption("show.signif.stars"), ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"), "\n\n", sep = "")
  cat("Pearson Residuals: \n")
  if (x$df.residual > 5) {
    x$pearson.resid = setNames(quantile(x$pearson.resid, na.rm = TRUE), c("Min", "1Q", "Median", "3Q", "Max"))
  }
  xx = zapsmall(x$pearson.resid, digits + 1L)
  print.default(xx, digits = digits, na.print = "", print.gap = 2L)
  if (length(x$aliased) == 0L) {
    cat("\nNo Coefficients\n")
  }
  else {
    df = if ("df" %in% names(x)) {
      x[["df"]]
    } else {NULL}
    if (!is.null(df) && (nsingular = df[3L] - df[1L])) {
      cat("\nCoefficients: (", nsingular, " not defined because of singularities)\n", sep = "")
    } else {cat("\nCoefficients:\n")}
    coefs = x$coefficients
    if (!is.null(aliased <- x$aliased) && any(aliased)) {
      cn = names(aliased)
      coefs = matrix(NA, length(aliased), 4L, dimnames = list(cn, colnames(coefs)))
      coefs[!aliased, ] = x$coefficients
    }
    printCoefmat(coefs, digits = digits, signif.stars = signif.stars, na.print = "NA", ...)
  }

  if(!is.null(x$deviance)){
    cat("\n", apply(cbind(paste(format(c("Null", "Residual"), justify = "right"), "deviance:"), format(unlist(x[c("null.deviance", "deviance")]), digits = max(5L, digits + 1L)), " on", format(unlist(x[c("df.null", "df.residual")])), " degrees of freedom\n"), 1L, paste, collapse = " "), sep = "")
    if (nzchar(mess <- naprint(x$na.action))) {cat("  (", mess, ")\n", sep = "")}
    cat("AIC:", format(x$aic, digits = max(4L, digits + 1L)), "\n")
    cat("UCMINF retries: ", x$retry$retries, ", iterations: ", x$iter[1], ", max gradient at convergence: ", format(x$iter[2], digits=max(5L, digits + 1L)), sep="")
  }

  correl = x$correlation
  if (!is.null(correl)) {
    p = NCOL(correl)
    if (p > 1) {
      cat("\nCorrelation of Coefficients:\n")
      if (is.logical(symbolic.cor) && symbolic.cor) {
        print(symnum(correl, abbr.colnames = NULL))
      } else {
        correl = format(round(correl, 2L), nsmall = 2L, digits = digits)
        correl[!lower.tri(correl)] = ""
        print(correl[-1, -p, drop = FALSE], quote = FALSE)
      }
    }
  }
  cat("\n")
  invisible(x)
}


## Calculate anova tables for spmle objects
#' @export
anova.spmle = function(object, ...) {
  ## If there's no loglikelihood, just print the object
  if(is.null(logLik(object))) {return(object)}

  ## If called on a single model, compare to the null model
  if(length(list(object, ...)) == 1L) {
    nullLik = -0.5*object$null.deviance
    attr(nullLik, "df") = 1
    null_model = list(nobs=object$nobs, logLik=nullLik, formula="Null Model")
    class(null_model) = "spmle"
    return(anova.spmle(null_model, object))
  }

  ## Collect all models specified as arguments
  objects = list(object, ...)
  nmodels = length(objects)
  ns = sapply(objects, nobs)
  if(any(ns != ns[1])) {stop("models were not all fitted to the same size of dataset")}

  ## Form ANOVA table
  rval = matrix(rep(NA, 5 * nmodels), ncol = 5)
  colnames(rval) = c("#Df", "LogLik", "Df", "Chisq", "Pr(>Chisq)")
  rownames(rval) = 1:nmodels
  logL = lapply(objects, logLik)

  ## Check that all are spmle models with loglikelihood
  if(!all(sapply(objects, inherits, "spmle"))) {stop("only models fit using spmle() may be compared")}
  if(any(sapply(logL, is.null))) {stop("models fit with spmleCombo() do not have an associated likelihood to be compared")}

  rval[, 1] = as.numeric(sapply(logL, function(x) attr(x, "df")))
  rval[, 2] = sapply(logL, as.numeric)
  rval[2:nmodels, 3] = rval[2:nmodels, 1] - rval[1:(nmodels - 1), 1]
  rval[2:nmodels, 4] = 2 * abs(rval[2:nmodels, 2] - rval[1:(nmodels - 1), 2])
  rval[, 5] = pchisq(rval[, 4], round(abs(rval[, 3])), lower.tail = FALSE)
  variables = lapply(objects, formula)
  title = "Likelihood ratio test\n"
  topnote = paste("Model ", format(1:nmodels), ": ", variables, sep = "", collapse = "\n")
  structure(as.data.frame(rval), heading = c(title, topnote), class = c("anova", "data.frame"))
}


## Plot spmle objects
#' @export
plot.spmle = function(x, ...) {
  resids = residuals(x)
  preds = fitted(x)
  frmula = as.character(as.expression(formula(x)))
  xlabel = bquote(paste("Predicted Values,  spmle: ", .(frmula), ",   ", pi[1]== .(x$pi1)))
  plot(x=preds, y=resids, xlab=xlabel, ylab="Pearson Residuals", main="Residuals vs Fitted")
  panel.smooth(preds, resids)
  abline(h = 0, lty = 3, col = "gray")
}


## Return covariance matrix
#' @export
vcov.spmle = function (object, ...) summary.spmle(object, ...)$cov.scaled


##  Calculate confidence intervals
#' @export
confint.spmle = function(object, parm, level=0.95, ...) {
  if(level>=1 & level<100) {level=level/100}
  stopifnot(0<level, level<1)

  pnames <- names(coef(object))
  if(missing(parm)) {
    parm <- seq_along(pnames)
  } else if(is.character(parm)) {
    parm <- match(parm, pnames, nomatch = 0L)
  }

  coefs = coef(object)[parm]
  SE = object$SE[parm]

  probs = c((1-level)/2, level+(1-level)/2)
  CI = cbind(qnorm(probs[1], mean=coefs, sd=SE), qnorm(probs[2], mean=coefs, sd=SE))
  colnames(CI) = paste(100*probs, "%")
  rownames(CI) = pnames[parm]
  return(CI)
}


## Extract model objects
#' @export
model.matrix.spmle = function(object, ...){
  model.matrix.lm(object$formula, object$model)
}

#' Predict method for spmle objects
#'
#' Obtains predictions from a fitted spmle object.
#'
#' \code{predict.spmle} produces predicted values, obtained by evaluating the
#' \code{spmle} function in the frame newdata (which defaults to \code{model.frame(object)}).
#' If the logical \code{se.fit} is \code{TRUE}, standard errors of the predictions
#' are calculated (only in the link scale). Setting \code{interval="confidence"} specifies computation of
#' confidence intervals at the specified level.
#'
#' @param object of class inheriting from "spmle"
#' @param newdata an optional list or data frame in which to look for variables
#'   to use when making predictions. If omitted, the fitted values are used.
#' @param se.fit a switch indicating if standard errors are required.
#' @param interval Type of interval calculation. Can be abbreviated. Prediction
#'   intervals are not meaningful for binary responses and are not allowed.
#' @param level confidence level.
#' @param type the type of prediction required. The default is \code{"link"} which
#'   uses the logit scale of the linear predictors, giving the log odds.
#'   The alternative \code{"response"} uses the probability scale, giving
#'   \code{Pr(D=1|G,E)}.
#' @param na.action function determining what should be done with missing values
#'   in newdata. The default is to predict NA.
#' @param ... further arguments passed to or from other methods.
#'
#' @return \code{predict.spmle} produces a vector of predictions or a matrix of
#'   predictions and bounds with column names fit, lwr, and upr if interval is set.
#'
#'   If \code{se.fit} is \code{TRUE}, a list with the following components is returned:
#'   \describe{
#'     \item{\code{fit}}{vector or matrix as above}
#'     \item{\code{se.fit}}{vector with standard error of predicted means, in the link scale}
#'     \item{\code{residual.scale}}{residual standard deviation}
#'     \item{\code{df}}{degrees of freedom for residual}
#'   }
#' @export
predict.spmle = function(object, newdata, se.fit=FALSE,
                         interval=c("none", "confidence"), level=0.95,
                         type=c("link", "response"), na.action=na.pass, ...) {
  if(level>=1 & level<100) {level=level/100}
  stopifnot(0<level, level<1)
  if(pmatch(interval[[1]], "prediction", nomatch=0)) {stop("prediction intervals for binary response are not meaningful")}

  tt <- terms(object)
  if(missing(newdata) || is.null(newdata)) {
    X <- model.matrix(object)
  } else {
    Terms <- delete.response(tt)
    m <- model.frame(Terms, newdata, na.action = na.action)
    if(!is.null(cl <- attr(Terms, "dataClasses"))) {.checkMFClasses(cl, m)}
    X <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
  }
  n <- length(object$residuals)
  p <- object$rank
  if (p < ncol(X) && !(missing(newdata) || is.null(newdata))) {warning("prediction from a rank-deficient fit may be misleading")}
  beta <- object$coefficients
  predictor <- drop(X %*% beta)

  interval <- match.arg(interval)
  type <- match.arg(type)
  if(se.fit && type=="response") {
    warning("standard errors can only be reported in the link scale")
    se.fit = FALSE
  }

  if (se.fit || interval != "none") {
    df <- object$df.residual
    res.var <-  sum(object$residuals^2)/df
    se = sqrt(diag(X %*% vcov(object) %*% t(X)))
  }
  if (interval != "none") {
    tfrac <- qt((1 - level)/2, df)
    hwid <- tfrac * se
    predictor <- cbind(predictor, predictor + hwid %o% c(1, -1))
    colnames(predictor) <- c("fit", "lwr", "upr")
  }
  if(type=="response") {predictor=plogis(predictor)}
  if(se.fit) {
    return(list(fit = predictor, se.fit = se, df = df, residual.scale = sqrt(res.var)))
  } else {return(predictor)}
}
