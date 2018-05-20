

# ox = optimx(par=Omega_start, fn=lik_fn, gr=grad_fn, hess=neghess, control=list(all.methods=TRUE, save.failures=TRUE, trace=0), D=D, G=G, E=E, pi1=pi1)

uc = ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)
nlr = nloptr(x0=Omega_start, eval_f=neglikgrad, opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel = 1e-4), D=D, G=G, E=E, pi1=pi1)


system.time(print(uc<-ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=list(invhessian.lt=solve(neghess(Omega=Omega_start, D=D, G=G, E=E, pi1=pi1))[lower.tri(diag(length(Omega_start)),diag=TRUE)]), D=D, G=G, E=E, pi1=pi1)))

system.time(print(uc<-ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)))
system.time(print(nlr<-nloptr(opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel=1e-4), x0=Omega_start, eval_f=neglikgrad, D=D, G=G, E=E, pi1=pi1)))
nlr$objective - uc$value
# 598.2329
# -1.130845 -0.006988747 0.2600060 0.4970546 -0.02545356 0.02253449 0.9529137 0.37286867 -0.03275638 -0.3224604 0.4583135 -0.2144142
NLOPT_LD_SLSQP, 3.8
NLOPT_LD_LBFGS, 2.3
NLOPT_LD_VAR1, 4.9
NLOPT_LD_VAR2, 4.9
NLOPT_LD_TNEWTON, 3.8
NLOPT_LD_TNEWTON_RESTART, 5
NLOPT_LD_TNEWTON_PRECOND, 3.6
NLOPT_LD_TNEWTON_PRECOND_RESTART, 3.2
NLOPT_LD_MMA, 7

################################################################################ ucminf control()
ucminf_con = list(trace=0, grtol=1e-06, xtol=1e-12, stepmax=1, maxeval=500,
                  grad="forward", gradstep=c(1e-06, 1e-08), invhessian.lt=NULL)
ucminf_con[names(ucminf_con) %in% names(con)] = con[names(con) %in% names(ucminf_con)]


################################################################################ z
> z=lm(1:10~sin(1:10))
> names(z)
[1] "coefficients"  "residuals"     "effects"       "rank"          "fitted.values" "assign"        "qr"            "df.residual"   "xlevels"       "call"
[11] "terms"         "model"

> str(z)
List of 12
$ coefficients : Named num [1:2] 5.6 -0.707
..- attr(*, "names")= chr [1:2] "(Intercept)" "sin(1:10)"
$ residuals    : Named num [1:10] -4.01 -2.96 -2.5 -2.13 -1.28 ...
..- attr(*, "names")= chr [1:10] "1" "2" "3" "4" ...
$ effects      : Named num [1:10] -17.39 1.55 -1.66 -1.93 -1.22 ...
..- attr(*, "names")= chr [1:10] "(Intercept)" "sin(1:10)" "" "" ...
$ rank         : int 2
$ fitted.values: Named num [1:10] 5.01 4.96 5.5 6.13 6.28 ...
..- attr(*, "names")= chr [1:10] "1" "2" "3" "4" ...
$ assign       : int [1:2] 0 1
$ qr           :List of 5
..$ qr   : num [1:10, 1:2] -3.162 0.316 0.316 0.316 0.316 ...
.. ..- attr(*, "dimnames")=List of 2
.. .. ..$ : chr [1:10] "1" "2" "3" "4" ...
.. .. ..$ : chr [1:2] "(Intercept)" "sin(1:10)"
.. ..- attr(*, "assign")= int [1:2] 0 1
..$ qraux: num [1:2] 1.32 1.27
..$ pivot: int [1:2] 1 2
..$ tol  : num 1e-07
..$ rank : int 2
..- attr(*, "class")= chr "qr"
$ df.residual  : int 8
$ xlevels      : Named list()
$ call         : language lm(formula = 1:10 ~ sin(1:10))
$ terms        :Classes 'terms', 'formula'  language 1:10 ~ sin(1:10)
.. ..- attr(*, "variables")= language list(1:10, sin(1:10))
.. ..- attr(*, "factors")= int [1:2, 1] 0 1
.. .. ..- attr(*, "dimnames")=List of 2
.. .. .. ..$ : chr [1:2] "1:10" "sin(1:10)"
.. .. .. ..$ : chr "sin(1:10)"
.. ..- attr(*, "term.labels")= chr "sin(1:10)"
.. ..- attr(*, "order")= int 1
.. ..- attr(*, "intercept")= int 1
.. ..- attr(*, "response")= int 1
.. ..- attr(*, ".Environment")=<environment: R_GlobalEnv>
  .. ..- attr(*, "predvars")= language list(1:10, sin(1:10))
.. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
.. .. ..- attr(*, "names")= chr [1:2] "1:10" "sin(1:10)"
$ model        :'data.frame':	10 obs. of  2 variables:
  ..$ 1:10     : int [1:10] 1 2 3 4 5 6 7 8 9 10
..$ sin(1:10): num [1:10] 0.841 0.909 0.141 -0.757 -0.959 ...
..- attr(*, "terms")=Classes 'terms', 'formula'  language 1:10 ~ sin(1:10)
.. .. ..- attr(*, "variables")= language list(1:10, sin(1:10))
.. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
.. .. .. ..- attr(*, "dimnames")=List of 2
.. .. .. .. ..$ : chr [1:2] "1:10" "sin(1:10)"
.. .. .. .. ..$ : chr "sin(1:10)"
.. .. ..- attr(*, "term.labels")= chr "sin(1:10)"
.. .. ..- attr(*, "order")= int 1
.. .. ..- attr(*, "intercept")= int 1
.. .. ..- attr(*, "response")= int 1
.. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv>
  .. .. ..- attr(*, "predvars")= language list(1:10, sin(1:10))
.. .. ..- attr(*, "dataClasses")= Named chr [1:2] "numeric" "numeric"
.. .. .. ..- attr(*, "names")= chr [1:2] "1:10" "sin(1:10)"
- attr(*, "class")= chr "lm"

################################################################################ lm
> lm
function (formula, data, subset, weights, na.action, method = "qr",
          model = TRUE, x = FALSE, y = FALSE, qr = TRUE, singular.ok = TRUE,
          contrasts = NULL, offset, ...)
{
  ret.x <- x
  ret.y <- y
  cl <- match.call()
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
               "offset"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE
  mf[[1L]] <- quote(stats::model.frame)
  mf <- eval(mf, parent.frame())
  if (method == "model.frame")
    return(mf)
  else if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
                     method), domain = NA)
  mt <- attr(mf, "terms")
  y <- model.response(mf, "numeric")
  w <- as.vector(model.weights(mf))
  if (!is.null(w) && !is.numeric(w))
    stop("'weights' must be a numeric vector")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(y))
      stop(gettextf("number of offsets is %d, should equal %d (number of observations)",
                    length(offset), NROW(y)), domain = NA)
  }
  if (is.empty.model(mt)) {
    x <- NULL
    z <- list(coefficients = if (is.matrix(y)) matrix(, 0,
                                                      3) else numeric(), residuals = y, fitted.values = 0 *
                y, weights = w, rank = 0L, df.residual = if (!is.null(w)) sum(w !=
                                                                                0) else if (is.matrix(y)) nrow(y) else length(y))
    if (!is.null(offset)) {
      z$fitted.values <- offset
      z$residuals <- y - offset
    }
  }
  else {
    x <- model.matrix(mt, mf, contrasts)
    z <- if (is.null(w))
      lm.fit(x, y, offset = offset, singular.ok = singular.ok,
             ...)
    else lm.wfit(x, y, w, offset = offset, singular.ok = singular.ok,
                 ...)
  }
  class(z) <- c(if (is.matrix(y)) "mlm", "lm")
  z$na.action <- attr(mf, "na.action")
  z$offset <- offset
  z$contrasts <- attr(x, "contrasts")
  z$xlevels <- .getXlevels(mt, mf)
  z$call <- cl
  z$terms <- mt
  if (model)
    z$model <- mf
  if (ret.x)
    z$x <- x
  if (ret.y)
    z$y <- y
  if (!qr)
    z$qr <- NULL
  z
}
<bytecode: 0x4f693a8>
  <environment: namespace:stats>

################################################################################ lm.fit
> lm.fit
function (x, y, offset = NULL, method = "qr", tol = 1e-07, singular.ok = TRUE,
          ...)
{
  if (is.null(n <- nrow(x)))
    stop("'x' must be a matrix")
  if (n == 0L)
    stop("0 (non-NA) cases")
  p <- ncol(x)
  if (p == 0L) {
    return(list(coefficients = numeric(), residuals = y,
                fitted.values = 0 * y, rank = 0, df.residual = length(y)))
  }
  ny <- NCOL(y)
  if (is.matrix(y) && ny == 1)
    y <- drop(y)
  if (!is.null(offset))
    y <- y - offset
  if (NROW(y) != n)
    stop("incompatible dimensions")
  if (method != "qr")
    warning(gettextf("method = '%s' is not supported. Using 'qr'",
                     method), domain = NA)
  chkDots(...)
  z <- .Call(C_Cdqrls, x, y, tol, FALSE)
  if (!singular.ok && z$rank < p)
    stop("singular fit encountered")
  coef <- z$coefficients
  pivot <- z$pivot
  r1 <- seq_len(z$rank)
  dn <- colnames(x)
  if (is.null(dn))
    dn <- paste0("x", 1L:p)
  nmeffects <- c(dn[pivot[r1]], rep.int("", n - z$rank))
  r2 <- if (z$rank < p)
    (z$rank + 1L):p
  else integer()
  if (is.matrix(y)) {
    coef[r2, ] <- NA
    if (z$pivoted)
      coef[pivot, ] <- coef
    dimnames(coef) <- list(dn, colnames(y))
    dimnames(z$effects) <- list(nmeffects, colnames(y))
  }
  else {
    coef[r2] <- NA
    if (z$pivoted)
      coef[pivot] <- coef
    names(coef) <- dn
    names(z$effects) <- nmeffects
  }
  z$coefficients <- coef
  r1 <- y - z$residuals
  if (!is.null(offset))
    r1 <- r1 + offset
  if (z$pivoted)
    colnames(z$qr) <- colnames(x)[z$pivot]
  qr <- z[c("qr", "qraux", "pivot", "tol", "rank")]
  c(z[c("coefficients", "residuals", "effects", "rank")], list(fitted.values = r1,
                                                               assign = attr(x, "assign"), qr = structure(qr, class = "qr"),
                                                               df.residual = n - z$rank))
}
<bytecode: 0x7158280>
  <environment: namespace:stats>

################################################################################ print.lm
> getAnywhere(print.lm)
A single object matching ‘print.lm’ was found
It was found in the following places
registered S3 method for print from namespace stats
namespace:stats
with value

function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  if (length(coef(x))) {
    cat("Coefficients:\n")
    print.default(format(coef(x), digits = digits), print.gap = 2L,
                  quote = FALSE)
  }
  else cat("No coefficients\n")
  cat("\n")
  invisible(x)
}
<bytecode: 0x5c38130>
  <environment: namespace:stats>

################################################################################ summary.lm
> getAnywhere(summary.lm)
A single object matching ‘summary.lm’ was found
It was found in the following places
package:stats
registered S3 method for summary from namespace stats
namespace:stats
with value

function (object, correlation = FALSE, symbolic.cor = FALSE,
          ...)
{
  z <- object
  p <- z$rank
  rdf <- z$df.residual
  if (p == 0) {
    r <- z$residuals
    n <- length(r)
    w <- z$weights
    if (is.null(w)) {
      rss <- sum(r^2)
    }
    else {
      rss <- sum(w * r^2)
      r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
    class(ans) <- "summary.lm"
    ans$aliased <- is.na(coef(object))
    ans$residuals <- r
    ans$df <- c(0L, n, length(ans$aliased))
    ans$coefficients <- matrix(NA, 0L, 4L)
    dimnames(ans$coefficients) <- list(NULL, c("Estimate",
                                               "Std. Error", "t value", "Pr(>|t|)"))
    ans$sigma <- sqrt(resvar)
    ans$r.squared <- ans$adj.r.squared <- 0
    return(ans)
  }
  if (is.null(z$terms))
    stop("invalid 'lm' object:  no 'terms' component")
  if (!inherits(object, "lm"))
    warning("calling summary.lm(<fake-lm-object>) ...")
  Qr <- qr.lm(object)
  n <- NROW(Qr$qr)
  if (is.na(z$df.residual) || n - p != z$df.residual)
    warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
  r <- z$residuals
  f <- z$fitted.values
  w <- z$weights
  if (is.null(w)) {
    mss <- if (attr(z$terms, "intercept"))
      sum((f - mean(f))^2)
    else sum(f^2)
    rss <- sum(r^2)
  }
  else {
    mss <- if (attr(z$terms, "intercept")) {
      m <- sum(w * f/sum(w))
      sum(w * (f - m)^2)
    }
    else sum(w * f^2)
    rss <- sum(w * r^2)
    r <- sqrt(w) * r
  }
  resvar <- rss/rdf
  if (is.finite(resvar) && resvar < (mean(f)^2 + var(f)) *
      1e-30)
    warning("essentially perfect fit: summary may be unreliable")
  p1 <- 1L:p
  R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
  se <- sqrt(diag(R) * resvar)
  est <- z$coefficients[Qr$pivot[p1]]
  tval <- est/se
  ans <- z[c("call", "terms", if (!is.null(z$weights)) "weights")]
  ans$residuals <- r
  ans$coefficients <- cbind(Estimate = est, `Std. Error` = se,
                            `t value` = tval, `Pr(>|t|)` = 2 * pt(abs(tval), rdf,
                                                                  lower.tail = FALSE))
  ans$aliased <- is.na(z$coefficients)
  ans$sigma <- sqrt(resvar)
  ans$df <- c(p, rdf, NCOL(Qr$qr))
  if (p != attr(z$terms, "intercept")) {
    df.int <- if (attr(z$terms, "intercept"))
      1L
    else 0L
    ans$r.squared <- mss/(mss + rss)
    ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
                                                       df.int)/rdf)
    ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
                        numdf = p - df.int, dendf = rdf)
  }
  else ans$r.squared <- ans$adj.r.squared <- 0
  ans$cov.unscaled <- R
  dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,
                                                             1)]
  if (correlation) {
    ans$correlation <- (R * resvar)/outer(se, se)
    dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
    ans$symbolic.cor <- symbolic.cor
  }
  if (!is.null(z$na.action))
    ans$na.action <- z$na.action
  class(ans) <- "summary.lm"
  ans
}
<bytecode: 0x723fba8>
  <environment: namespace:stats>
