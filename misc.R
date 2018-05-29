
set.seed(1)
dat = simulate_complex(ncase=200,
                       ncontrol=200,
                       beta0=-4.14,
                       betaG_SNP=c(log(1.2), log(1.2), 0, log(1.2), 0),
                       betaE_bin=log(1.5),
                       betaGE_SNP_bin=c(log(1.3), 0, 0, log(1.3), 0),
                       MAF=c(0.1, 0.3, 0.3, 0.3, 0.1),
                       SNP_cor=0.7,
                       E_bin_freq=0.5)

z=spmle(D,G,E,.03,dat)
z
summary(z)

zz=combo_spmle(D,G,E,.03,dat,20)
zz
summary(zz)
head(predict(zz, interval="conf"))

zzz=combo_spmle(D,G,E,pi1=.03,data=dat,nboot=20, ncores=2, control=list(trace=1))
zzz
summary(zzz)

Z=combo_spmle(D,G,E,pi1=.03,data=dat,nboot=100, ncores=4)
Z
summary(Z)

################################################################################ Adding data to spmle
## Calculate asymptotic SE for SPMLE
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

  ## Save model frame
  model = model.frame(formula=formula, data=data)

  ## If starting values weren't provided, use logistic regression estimates
  if(missing(startvals)) {
    Omega_start = coef(glm(formula, family=binomial(link='logit'), data=data))
  } else {
    Omega_start = startvals
  }

  ## If user-provided startval lacked names, add them
  if(is.null(names(Omega_start))) {names(Omega_start) = colnames(model.matrix(formula, model[1,]))}

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
    swap_order = c(1,(2+nE):(1+nG+nE),2:(1+nE), rep(seq(from=(2+nG+nE), by=nE, length.out=nG), times=nE) + rep(0:(nE-1), each=nG))
    reverse_swap_order = c(1,(2+nG):(1+nG+nE),2:(1+nG), rep(seq(from=(2+nG+nE), by=nG, length.out=nE), times=nG) + rep(0:(nG-1), each=nE))
    temp = G
    G = E
    E = temp
  } else {
    swap_order = reverse_swap_order = seq_along(Omega_start)
  }

  ## Precondition with hessian if requested
  if(con$use_hess==TRUE) {con$invhessian.lt = solve(neghess(Omega=Omega_start, D=D, G=G, E=E, pi1=pi1))[lower.tri(diag(length(Omega_start)),diag=TRUE)]}

  ## Calculate SPMLE for a given disease rate
  spmle_max = maximize_spmle(Omega_start=Omega_start[reverse_swap_order], D=D, G=G, E=E, pi1=pi1, control=con)

  ## Asymptotic SE for known disease rate
  hess_zeta = hesszeta(Omega=spmle_max$par, D=D, G=G, E=E, pi1=pi1)
  Sigma = ((ncontrol-1)*cov(hess_zeta$zeta0) + (ncase-1)*cov(hess_zeta$zeta1))/n  # (ncontrol-1) to correct denominator because cov uses the "sample covariance matrix" denominator of n-1
  H_inv = -n*solve(hess_zeta$hessian)    # = (-hessian/n)^-1 = (Gamma_1 - Gamma_2)^-1
  Lambda = H_inv %*% Sigma %*% t(H_inv)  # covar matrix of sqrt(n) * OmegaHat
  SE_asy = sqrt(diag(Lambda)/n)

  ## Compile results into a list.  If we're swapping G & E, change back now
  spmle_est = list(coefficients = spmle_max$par[swap_order],
                   SE    = SE_asy[swap_order],
                   cov   = Lambda[swap_order,swap_order]/n,
                   H_inv = H_inv[swap_order,swap_order],
                   Sigma = Sigma[swap_order,swap_order],
                   zeta0 = hess_zeta$zeta0[, swap_order],
                   zeta1 = hess_zeta$zeta1[, swap_order],
                   ucminf = spmle_max,
                   call  = cl,
                   formula = formula,
                   data = data,
                   model = model)

  ## Return an object of class spmle, which inherits from glm and lm
  class(spmle_est) = c("spmle")#, "glm", "lm")
  return(spmle_est)
}


################################################################################ How to handle data frames (from CatCoxMLE)
## Get variables from data frame, if necessary
if(class(data)=="matrix") data = as.data.frame(data)
measuredTime = as.double(with(data=data, eval(measuredTimeName)))
failIndicator = as.integer(with(data=data, eval(failIndicatorName)))
probCat = as.matrix(with(data=data, eval(probCatName)))
covars = as.matrix(with(data=data, eval(covarsName)))

## Save their column names
probCatNames = colnames(probCat, do.NULL=FALSE, prefix=as.character(probCatName))
covarsNames = colnames(covars, do.NULL=FALSE, prefix=as.character(covarsName))

unconstrainedPar = setNames(object=unconstrainedMLE$par, nm=c(probCatNames[-1], covarsNames))


################################################################################ Checking what to include in the value of spmle return
set.seed(11)
x=rnorm(9)
# xx=rnorm(8)
y=x+rnorm(9, sd=3)>0
z=glm(y~x, family=binomial)
head(model.matrix(z))
head(model.frame(z))

summary(z$fitted.values)
summary(z$residuals)
sum(-2*log(abs(1-y-z$fitted.values)))

names(z)

sp = spmle(D=dat$D, G=dat$G, E=dat$E, pi1=0.3)
sp = spmle(D=D, G=G, E=E, pi1=pi1)
mm = model.matrix(~G*E)
linear_predictors = mm %*% sp$par

G[1,] = E[1,] = 10

set.seed(11)
x=rexp(9, .01)
y=rpois(n=9, x)^2
xydat = data.frame(xpart=x, ypart=y, silly=0)
z=glm(ypart~xpart, family=quasipoisson, data=xydat)
z=glm(ypart~exp(xpart), family=quasipoisson, data=xydat)
summary(z)
################################################################################ glm methods
> methods(class=glm)
[1] add1           anova          coerce         confint        cooks.distance deviance       drop1          effects        extractAIC     family         formula
[12] influence      initialize     logLik         model.frame    nobs           predict        print          residuals      rstandard      rstudent       show
[23] slotsFromS3    summary        vcov           weights


################################################################################ glm names
> names(z)
[1] "coefficients"      "residuals"         "fitted.values"     "effects"           "R"                 "rank"              "qr"                "family"
[9] "linear.predictors" "deviance"          "aic"               "null.deviance"     "iter"              "weights"           "prior.weights"     "df.residual"
[17] "df.null"           "y"                 "converged"         "boundary"          "model"             "call"              "formula"           "terms"
[25] "data"              "offset"            "control"           "method"            "contrasts"         "xlevels"

################################################################################ simulate & attach data
Rcpp::sourceCpp("~/Research/Tianying/caseControlGE/src/lik_grad_hess_zeta_02.cpp")
source("~/Research/Tianying/caseControlGE/R/simulate_complex.R")
set.seed(1)
dat = simulate_complex(ncase=500,
                       ncontrol=500,
                       beta0=-3.95,
                       betaG_SNP=c(log(1.2), log(1.2), 0),
                       betaE_bin=c(log(1.5), 0),
                       betaGE_SNP_bin=c(log(1.3), 0, 0, log(1.3), 0, 0),
                       MAF=c(0.1, 0.3, 0.3),
                       SNP_cor=0.7,
                       E_bin_freq=0.5)
D=dat$D ; G=dat$G ; E=dat$E ; pi1 = 0.03
################################################################################ loop to rerun optimization if convergence is not met

if(!is.finite(SPMLE$info[1]) | SPMLE$info[1]>con$max_grad_tol) {
  ## Scale starting values by 1/SD
  scale_vec = 1/c(1, apply(model.matrix(~0+G*E), 2, sd))

  ## Loop through different starting values
  for(j in seq_len(con$num_retries)) {
    if(con$trace>-1) cat("UCMINF retry", j, "of", con$num_retries, "\n")

    ## If failure happened when preconditioning with the hessian, try without (same start values the first time)
    if(!is.null(ucminf_con$invhessian.lt)){
      ucminf_con$invhessian.lt = NULL
      SPMLE = ucminf::ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)
    } else {
      SPMLE = ucminf::ucminf(par=rnorm(n=length(Omega_start), sd=scale_vec), fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)
    }

    ## Break out of the loop once we have convergence
    if(is.finite(SPMLE$info[1]) & SPMLE$info[1]<con$max_grad_tol) break
  }

  ## Test again: if still no convergance, throw an error
  if(!is.finite(SPMLE$info[1]) | SPMLE$info[1]>con$max_grad_tol) stop("UCMINF failed to converge")
}


################################################################################ test optimization algorithms

# ox = optimx(par=Omega_start, fn=lik_fn, gr=grad_fn, hess=neghess, control=list(all.methods=TRUE, save.failures=TRUE, trace=0), D=D, G=G, E=E, pi1=pi1)

uc = ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)
nlr = nloptr(x0=Omega_start, eval_f=neglikgrad, opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel = 1e-4), D=D, G=G, E=E, pi1=pi1)

system.time(print(uc<-ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=list(invhessian.lt=solve(neghess(Omega=Omega_start, D=D, G=G, E=E, pi1=pi1))[lower.tri(diag(length(Omega_start)),diag=TRUE)]), D=D, G=G, E=E, pi1=pi1)))

system.time(print(uc<-ucminf(par=Omega_start, fn=lik_fn, gr=grad_fn, control=ucminf_con, D=D, G=G, E=E, pi1=pi1)))
system.time(print(nlr<-nloptr(opts=list(algorithm="NLOPT_LD_LBFGS", xtol_rel=1e-4), x0=Omega_start, eval_f=neglikgrad, D=D, G=G, E=E, pi1=pi1)))
nlr$objective - uc$value
# 598.2329
# -1.130845 -0.006988747 0.2600060 0.4970546 -0.02545356 0.02253449 0.9529137 0.37286867 -0.03275638 -0.3224604 0.4583135 -0.2144142
# NLOPT_LD_SLSQP, 3.8
# NLOPT_LD_LBFGS, 2.3
# NLOPT_LD_VAR1, 4.9
# NLOPT_LD_VAR2, 4.9
# NLOPT_LD_TNEWTON, 3.8
# NLOPT_LD_TNEWTON_RESTART, 5
# NLOPT_LD_TNEWTON_PRECOND, 3.6
# NLOPT_LD_TNEWTON_PRECOND_RESTART, 3.2
# NLOPT_LD_MMA, 7 (didn't converge)

################################################################################ ucminf control()

ucminf_con_names = c("trace", "grtol", "xtol", "stepmax", "maxeval", "grad", "gradstep", "invhessian.lt")
ucminf_con = con[names(con) %in% ucminf_con_names]


ucminf_con = list(trace=0, grtol=1e-06, xtol=1e-12, stepmax=1, maxeval=500,
                  grad="forward", gradstep=c(1e-06, 1e-08), invhessian.lt=NULL)
ucminf_con[names(ucminf_con) %in% names(con)] = con[names(con) %in% names(ucminf_con)]


################################################################################ z glm
set.seed(1)
x=rnorm(200)
y=x+rnorm(200, .1)>0
z=glm(y~x, family=binomial)

> str(z)
List of 30
 $ coefficients     : Named num [1:2] 0.105 1.407
  ..- attr(*, "names")= chr [1:2] "(Intercept)" "x"
 $ residuals        : Named num [1:200] -1.46 1.7 3.92 1.1 -2.76 ...
  ..- attr(*, "names")= chr [1:200] "1" "2" "3" "4" ...
 $ fitted.values    : Named num [1:200] 0.315 0.59 0.255 0.913 0.638 ...
  ..- attr(*, "names")= chr [1:200] "1" "2" "3" "4" ...
 $ effects          : Named num [1:200] -0.322 -6.108 1.818 0.245 -1.318 ...
  ..- attr(*, "names")= chr [1:200] "(Intercept)" "x" "" "" ...
 $ R                : num [1:2, 1:2] -6.128 0 0.227 -4.342
  ..- attr(*, "dimnames")=List of 2
  .. ..$ : chr [1:2] "(Intercept)" "x"
  .. ..$ : chr [1:2] "(Intercept)" "x"
 $ rank             : int 2
 $ qr               :List of 5
  ..$ qr   : num [1:200, 1:2] -6.1278 0.0803 0.0711 0.046 0.0784 ...
  .. ..- attr(*, "dimnames")=List of 2
  .. .. ..$ : chr [1:200] "1" "2" "3" "4" ...
  .. .. ..$ : chr [1:2] "(Intercept)" "x"
  ..$ rank : int 2
  ..$ qraux: num [1:2] 1.08 1.03
  ..$ pivot: int [1:2] 1 2
  ..$ tol  : num 1e-11
  ..- attr(*, "class")= chr "qr"
 $ family           :List of 12
  ..$ family    : chr "binomial"
  ..$ link      : chr "logit"
  ..$ linkfun   :function (mu)
  ..$ linkinv   :function (eta)
  ..$ variance  :function (mu)
  ..$ dev.resids:function (y, mu, wt)
  ..$ aic       :function (y, n, mu, wt, dev)
  ..$ mu.eta    :function (eta)
  ..$ initialize:  expression({  if (NCOL(y) == 1) {  if (is.factor(y))  y <- y != levels(y)[1L]  n <- rep.int(1, nobs)  y[weights =| __truncated__
  ..$ validmu   :function (mu)
  ..$ valideta  :function (eta)
  ..$ simulate  :function (object, nsim)
  ..- attr(*, "class")= chr "family"
 $ linear.predictors: Named num [1:200] -0.777 0.363 -1.071 2.349 0.568 ...
  ..- attr(*, "names")= chr [1:200] "1" "2" "3" "4" ...
 $ deviance         : num 222
 $ aic              : num 226
 $ null.deviance    : num 277
 $ iter             : int 4
 $ weights          : Named num [1:200] 0.2158 0.2419 0.1901 0.0796 0.2309 ...
  ..- attr(*, "names")= chr [1:200] "1" "2" "3" "4" ...
 $ prior.weights    : Named num [1:200] 1 1 1 1 1 1 1 1 1 1 ...
  ..- attr(*, "names")= chr [1:200] "1" "2" "3" "4" ...
 $ df.residual      : int 198
 $ df.null          : int 199
 $ y                : Named num [1:200] 0 1 1 1 0 1 1 1 1 1 ...
  ..- attr(*, "names")= chr [1:200] "1" "2" "3" "4" ...
 $ converged        : logi TRUE
 $ boundary         : logi FALSE
 $ model            :'data.frame':	200 obs. of  2 variables:
  ..$ y: logi [1:200] FALSE TRUE TRUE TRUE FALSE TRUE ...
  ..$ x: num [1:200] -0.626 0.184 -0.836 1.595 0.33 ...
  ..- attr(*, "terms")=Classes 'terms', 'formula'  language y ~ x
  .. .. ..- attr(*, "variables")= language list(y, x)
  .. .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. .. ..$ : chr [1:2] "y" "x"
  .. .. .. .. ..$ : chr "x"
  .. .. ..- attr(*, "term.labels")= chr "x"
  .. .. ..- attr(*, "order")= int 1
  .. .. ..- attr(*, "intercept")= int 1
  .. .. ..- attr(*, "response")= int 1
  .. .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv>
  .. .. ..- attr(*, "predvars")= language list(y, x)
  .. .. ..- attr(*, "dataClasses")= Named chr [1:2] "logical" "numeric"
  .. .. .. ..- attr(*, "names")= chr [1:2] "y" "x"
 $ call             : language glm(formula = y ~ x, family = binomial)
 $ formula          :Class 'formula'  language y ~ x
  .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv>
 $ terms            :Classes 'terms', 'formula'  language y ~ x
  .. ..- attr(*, "variables")= language list(y, x)
  .. ..- attr(*, "factors")= int [1:2, 1] 0 1
  .. .. ..- attr(*, "dimnames")=List of 2
  .. .. .. ..$ : chr [1:2] "y" "x"
  .. .. .. ..$ : chr "x"
  .. ..- attr(*, "term.labels")= chr "x"
  .. ..- attr(*, "order")= int 1
  .. ..- attr(*, "intercept")= int 1
  .. ..- attr(*, "response")= int 1
  .. ..- attr(*, ".Environment")=<environment: R_GlobalEnv>
  .. ..- attr(*, "predvars")= language list(y, x)
  .. ..- attr(*, "dataClasses")= Named chr [1:2] "logical" "numeric"
  .. .. ..- attr(*, "names")= chr [1:2] "y" "x"
 $ data             :<environment: R_GlobalEnv>
 $ offset           : NULL
 $ control          :List of 3
  ..$ epsilon: num 1e-08
  ..$ maxit  : num 25
  ..$ trace  : logi FALSE
 $ method           : chr "glm.fit"
 $ contrasts        : NULL
 $ xlevels          : Named list()
 - attr(*, "class")= chr [1:2] "glm" "lm"

################################################################################ glm
> glm
function (formula, family = gaussian, data, weights, subset,
    na.action, start = NULL, etastart, mustart, offset, control = list(...),
    model = TRUE, method = "glm.fit", x = FALSE, y = TRUE, contrasts = NULL,
    ...)
{
    call <- match.call()
    if (is.character(family))
        family <- get(family, mode = "function", envir = parent.frame())
    if (is.function(family))
        family <- family()
    if (is.null(family$family)) {
        print(family)
        stop("'family' not recognized")
    }
    if (missing(data))
        data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- quote(stats::model.frame)
    mf <- eval(mf, parent.frame())
    if (identical(method, "model.frame"))
        return(mf)
    if (!is.character(method) && !is.function(method))
        stop("invalid 'method' argument")
    if (identical(method, "glm.fit"))
        control <- do.call("glm.control", control)
    mt <- attr(mf, "terms")
    Y <- model.response(mf, "any")
    if (length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if (!is.null(nm))
            names(Y) <- nm
    }
    X <- if (!is.empty.model(mt))
        model.matrix(mt, mf, contrasts)
    else matrix(, NROW(Y), 0L)
    weights <- as.vector(model.weights(mf))
    if (!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if (!is.null(weights) && any(weights < 0))
        stop("negative weights not allowed")
    offset <- as.vector(model.offset(mf))
    if (!is.null(offset)) {
        if (length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)",
                length(offset), NROW(Y)), domain = NA)
    }
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
    fit <- eval(call(if (is.function(method)) "method" else method,
        x = X, y = Y, weights = weights, start = start, etastart = etastart,
        mustart = mustart, offset = offset, family = family,
        control = control, intercept = attr(mt, "intercept") >
            0L))
    if (length(offset) && attr(mt, "intercept") > 0L) {
        fit2 <- eval(call(if (is.function(method)) "method" else method,
            x = X[, "(Intercept)", drop = FALSE], y = Y, weights = weights,
            offset = offset, family = family, control = control,
            intercept = TRUE))
        if (!fit2$converged)
            warning("fitting to calculate the null deviance did not converge -- increase 'maxit'?")
        fit$null.deviance <- fit2$deviance
    }
    if (model)
        fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if (x)
        fit$x <- X
    if (!y)
        fit$y <- NULL
    fit <- c(fit, list(call = call, formula = formula, terms = mt,
        data = data, offset = offset, control = control, method = method,
        contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt,
            mf)))
    class(fit) <- c(fit$class, c("glm", "lm"))
    fit
}
<bytecode: 0x473ca70>
<environment: namespace:stats>

################################################################################ glm.fit
> glm.fit
function (x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
    mustart = NULL, offset = rep(0, nobs), family = gaussian(),
    control = list(), intercept = TRUE)
{
    control <- do.call("glm.control", control)
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    linkinv <- family$linkinv
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object",
            call. = FALSE)
    dev.resids <- family$dev.resids
    aic <- family$aic
    mu.eta <- family$mu.eta
    unless.null <- function(x, if.null) if (is.null(x))
        if.null
    else x
    valideta <- unless.null(family$valideta, function(eta) TRUE)
    validmu <- unless.null(family$validmu, function(mu) TRUE)
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model",
                call. = FALSE)
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model", call. = FALSE)
        dev <- sum(dev.resids(y, mu, weights))
        w <- sqrt((weights * mu.eta(eta)^2)/variance(mu))
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep_len(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric()
        iter <- 0L
    }
    else {
        coefold <- NULL
        eta <- if (!is.null(etastart))
            etastart
        else if (!is.null(start))
            if (length(start) != nvars)
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                  nvars, paste(deparse(xnames), collapse = ", ")),
                  domain = NA)
            else {
                coefold <- start
                offset + as.vector(if (NCOL(x) == 1L)
                  x * start
                else x %*% start)
            }
        else family$linkfun(mustart)
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some",
                call. = FALSE)
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        for (iter in 1L:control$maxit) {
            good <- weights > 0
            varmu <- variance(mu)[good]
            if (anyNA(varmu))
                stop("NAs in V(mu)")
            if (any(varmu == 0))
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            if (any(is.na(mu.eta.val[good])))
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning(gettextf("no observations informative at iteration %d",
                  iter), domain = NA)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/variance(mu)[good])
            fit <- .Call(C_Cdqrls, x[good, , drop = FALSE] *
                w, z * w, min(1e-07, control$epsilon/1000), check = FALSE)
            if (any(!is.finite(fit$coefficients))) {
                conv <- FALSE
                warning(gettextf("non-finite coefficients at iteration %d",
                  iter), domain = NA)
                break
            }
            if (nobs < fit$rank)
                stop(sprintf(ngettext(nobs, "X matrix has rank %d, but only %d observation",
                  "X matrix has rank %d, but only %d observations"),
                  fit$rank, nobs), domain = NA)
            start[fit$pivot] <- fit$coefficients
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)
                cat("Deviance = ", dev, " Iterations - ", iter,
                  "\n", sep = "")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated due to divergence",
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                  cat("Step halved: new deviance = ", dev, "\n",
                    sep = "")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated: out of bounds",
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size",
                      call. = FALSE)
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace)
                  cat("Step halved: new deviance = ", dev, "\n",
                    sep = "")
            }
            if (abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {
                devold <- dev
                coef <- coefold <- start
            }
        }
        if (!conv)
            warning("glm.fit: algorithm did not converge", call. = FALSE)
        if (boundary)
            warning("glm.fit: algorithm stopped at boundary value",
                call. = FALSE)
        eps <- 10 * .Machine$double.eps
        if (family$family == "binomial") {
            if (any(mu > 1 - eps) || any(mu < eps))
                warning("glm.fit: fitted probabilities numerically 0 or 1 occurred",
                  call. = FALSE)
        }
        if (family$family == "poisson") {
            if (any(mu < eps))
                warning("glm.fit: fitted rates numerically 0 occurred",
                  call. = FALSE)
        }
        if (fit$rank < nvars)
            coef[fit$pivot][seq.int(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- (y - mu)/mu.eta(eta)
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1L:nr, 1L:nvars] <- fit$qr[1L:nr, 1L:nvars]
        }
        else Rmat <- fit$qr[1L:nvars, 1L:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    }
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    if (!EMPTY)
        names(fit$effects) <- c(xxnames[seq_len(fit$rank)], rep.int("",
            sum(good) - fit$rank))
    wtdmu <- if (intercept)
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY)
        0
    else fit$rank
    resdf <- n.ok - rank
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    list(coefficients = coef, residuals = residuals, fitted.values = mu,
        effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
        rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank",
            "qraux", "pivot", "tol")], class = "qr"), family = family,
        linear.predictors = eta, deviance = dev, aic = aic.model,
        null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
        df.residual = resdf, df.null = nulldf, y = y, converged = conv,
        boundary = boundary)
}
<bytecode: 0x4fe5f40>
<environment: namespace:stats>


################################################################################ glm methods
> methods(class=glm)
[1] add1           anova          coerce         confint        cooks.distance deviance       drop1          effects        extractAIC     family         formula
[12] influence      initialize     logLik         model.frame    nobs           predict        print          residuals      rstandard      rstudent       show
[23] slotsFromS3    summary        vcov           weights


################################################################################ summary.glm
> getAnywhere(summary.glm)
A single object matching ‘summary.glm’ was found
It was found in the following places
  package:stats
  registered S3 method for summary from namespace stats
  namespace:stats
with value

function (object, dispersion = NULL, correlation = FALSE, symbolic.cor = FALSE,
    ...)
{
    est.disp <- FALSE
    df.r <- object$df.residual
    if (is.null(dispersion)) {dispersion <- if (object$family$family %in% c("poisson", "binomial")) 1
        } else if (df.r > 0) {
            est.disp <- TRUE
            if (any(object$weights == 0)) warning("observations with zero weight not used for calculating dispersion")
            sum((object$weights * object$residuals^2)[object$weights >0])/df.r
        } else {
            est.disp <- TRUE
            NaN
        }
    aliased <- is.na(coef(object))
    p <- object$rank
    if (p > 0) {
        p1 <- 1L:p
        Qr <- stats:::qr.lm(object)
        coef.p <- object$coefficients[Qr$pivot[p1]]
        covmat.unscaled <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
        dimnames(covmat.unscaled) <- list(names(coef.p), names(coef.p))
        covmat <- dispersion * covmat.unscaled
        var.cf <- diag(covmat)
        s.err <- sqrt(var.cf)
        tvalue <- coef.p/s.err
        dn <- c("Estimate", "Std. Error")
        if (!est.disp) {
            pvalue <- 2 * pnorm(-abs(tvalue))
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                "z value", "Pr(>|z|)"))
        }
        else if (df.r > 0) {
            pvalue <- 2 * pt(-abs(tvalue), df.r)
            coef.table <- cbind(coef.p, s.err, tvalue, pvalue)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                "t value", "Pr(>|t|)"))
        }
        else {
            coef.table <- cbind(coef.p, NaN, NaN, NaN)
            dimnames(coef.table) <- list(names(coef.p), c(dn,
                "t value", "Pr(>|t|)"))
        }
        df.f <- NCOL(Qr$qr)
    }
    else {
        coef.table <- matrix(, 0L, 4L)
        dimnames(coef.table) <- list(NULL, c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)"))
        covmat.unscaled <- covmat <- matrix(, 0L, 0L)
        df.f <- length(aliased)
    }
    keep <- match(c("call", "terms", "family", "deviance", "aic",
        "contrasts", "df.residual", "null.deviance", "df.null",
        "iter", "na.action"), names(object), 0L)
    ans <- c(object[keep], list(deviance.resid = residuals(object,
        type = "deviance"), coefficients = coef.table, aliased = aliased,
        dispersion = dispersion, df = c(object$rank, df.r, df.f),
        cov.unscaled = covmat.unscaled, cov.scaled = covmat))
    if (correlation && p > 0) {
        dd <- sqrt(diag(covmat.unscaled))
        ans$correlation <- covmat.unscaled/outer(dd, dd)
        ans$symbolic.cor <- symbolic.cor
    }
    class(ans) <- "summary.glm"
    return(ans)
}
<bytecode: 0x85bff58>
<environment: namespace:stats>


################################################################################ print.glm
> getAnywhere(print.glm)
A single object matching ‘print.glm’ was found
It was found in the following places
  registered S3 method for print from namespace stats
  namespace:stats
with value

function (x, digits = max(3L, getOption("digits") - 3L), ...)
{
    cat("\nCall:  ", paste(deparse(x$call), sep = "\n", collapse = "\n"),
        "\n\n", sep = "")
    if (length(coef(x))) {
        cat("Coefficients")
        if (is.character(co <- x$contrasts))
            cat("  [contrasts: ", apply(cbind(names(co), co),
                1L, paste, collapse = "="), "]")
        cat(":\n")
        print.default(format(x$coefficients, digits = digits),
            print.gap = 2, quote = FALSE)
    }
    else cat("No coefficients\n\n")
    cat("\nDegrees of Freedom:", x$df.null, "Total (i.e. Null); ",
        x$df.residual, "Residual\n")
    if (nzchar(mess <- naprint(x$na.action)))
        cat("  (", mess, ")\n", sep = "")
    cat("Null Deviance:\t   ", format(signif(x$null.deviance,
        digits)), "\nResidual Deviance:", format(signif(x$deviance,
        digits)), "\tAIC:", format(signif(x$aic, digits)))
    cat("\n")
    invisible(x)
}
<bytecode: 0x85fc428>
<environment: namespace:stats>


################################################################################ anova.glm
> getAnywhere(anova.glm)
A single object matching ‘anova.glm’ was found
It was found in the following places
  registered S3 method for anova from namespace stats
  namespace:stats
with value

function (object, ..., dispersion = NULL, test = NULL)
{
    dotargs <- list(...)
    named <- if (is.null(names(dotargs)))
        rep_len(FALSE, length(dotargs))
    else (names(dotargs) != "")
    if (any(named))
        warning("the following arguments to 'anova.glm' are invalid and dropped: ",
            paste(deparse(dotargs[named]), collapse = ", "))
    dotargs <- dotargs[!named]
    is.glm <- vapply(dotargs, function(x) inherits(x, "glm"),
        NA)
    dotargs <- dotargs[is.glm]
    if (length(dotargs))
        return(anova.glmlist(c(list(object), dotargs), dispersion = dispersion,
            test = test))
    doscore <- !is.null(test) && test == "Rao"
    varlist <- attr(object$terms, "variables")
    x <- if (n <- match("x", names(object), 0L))
        object[[n]]
    else model.matrix(object)
    varseq <- attr(x, "assign")
    nvars <- max(0, varseq)
    resdev <- resdf <- NULL
    if (doscore) {
        score <- numeric(nvars)
        method <- object$method
        y <- object$y
        fit <- eval(call(if (is.function(method)) "method" else method,
            x = x[, varseq == 0, drop = FALSE], y = y, weights = object$prior.weights,
            start = object$start, offset = object$offset, family = object$family,
            control = object$control))
        r <- fit$residuals
        w <- fit$weights
    }
    if (nvars > 1 || doscore) {
        method <- object$method
        y <- object$y
        if (is.null(y)) {
            mu.eta <- object$family$mu.eta
            eta <- object$linear.predictors
            y <- object$fitted.values + object$residuals * mu.eta(eta)
        }
        for (i in seq_len(nvars - 1L)) {
            fit <- eval(call(if (is.function(method)) "method" else method,
                x = x[, varseq <= i, drop = FALSE], y = y, weights = object$prior.weights,
                start = object$start, offset = object$offset,
                family = object$family, control = object$control))
            if (doscore) {
                zz <- eval(call(if (is.function(method)) "method" else method,
                  x = x[, varseq <= i, drop = FALSE], y = r,
                  weights = w))
                score[i] <- zz$null.deviance - zz$deviance
                r <- fit$residuals
                w <- fit$weights
            }
            resdev <- c(resdev, fit$deviance)
            resdf <- c(resdf, fit$df.residual)
        }
        if (doscore) {
            zz <- eval(call(if (is.function(method)) "method" else method,
                x = x, y = r, weights = w))
            score[nvars] <- zz$null.deviance - zz$deviance
        }
    }
    resdf <- c(object$df.null, resdf, object$df.residual)
    resdev <- c(object$null.deviance, resdev, object$deviance)
    table <- data.frame(c(NA, -diff(resdf)), c(NA, pmax(0, -diff(resdev))),
        resdf, resdev)
    tl <- attr(object$terms, "term.labels")
    if (length(tl) == 0L)
        table <- table[1, , drop = FALSE]
    dimnames(table) <- list(c("NULL", tl), c("Df", "Deviance",
        "Resid. Df", "Resid. Dev"))
    if (doscore)
        table <- cbind(table, Rao = c(NA, score))
    title <- paste0("Analysis of Deviance Table", "\n\nModel: ",
        object$family$family, ", link: ", object$family$link,
        "\n\nResponse: ", as.character(varlist[-1L])[1L], "\n\nTerms added sequentially (first to last)\n\n")
    df.dispersion <- Inf
    if (is.null(dispersion)) {
        dispersion <- summary(object, dispersion = dispersion)$dispersion
        df.dispersion <- if (dispersion == 1)
            Inf
        else object$df.residual
    }
    if (!is.null(test)) {
        if (test == "F" && df.dispersion == Inf) {
            fam <- object$family$family
            if (fam == "binomial" || fam == "poisson")
                warning(gettextf("using F test with a '%s' family is inappropriate",
                  fam), domain = NA)
            else warning("using F test with a fixed dispersion is inappropriate")
        }
        table <- stat.anova(table = table, test = test, scale = dispersion,
            df.scale = df.dispersion, n = NROW(x))
    }
    structure(table, heading = title, class = c("anova", "data.frame"))
}
<bytecode: 0x8624928>
<environment: namespace:stats>


################################################################################ predict.glm
> getAnywhere(predict.glm)
A single object matching ‘predict.glm’ was found
It was found in the following places
  package:stats
  registered S3 method for predict from namespace stats
  namespace:stats
with value

predict.glm = function(object, newdata = NULL, type = c("link", "response",
    "terms"), se.fit = FALSE, dispersion = NULL, terms = NULL,
    na.action = na.pass, ...)
{
    type <- match.arg(type)
    na.act <- object$na.action
    object$na.action <- NULL
    if (!se.fit) {
        if (missing(newdata)) {
            pred <- switch(type, link = object$linear.predictors,
                response = object$fitted.values, terms = predict.lm(object,
                  se.fit = se.fit, scale = 1, type = "terms",
                  terms = terms))
            if (!is.null(na.act))
                pred <- napredict(na.act, pred)
        }
        else {
            pred <- predict.lm(object, newdata, se.fit, scale = 1,
                type = ifelse(type == "link", "response", type),
                terms = terms, na.action = na.action)
            switch(type, response = {
                pred <- family(object)$linkinv(pred)
            }, link = , terms = )
        }
    }
    else {
        if (inherits(object, "survreg"))
            dispersion <- 1
        if (is.null(dispersion) || dispersion == 0)
            dispersion <- summary(object, dispersion = dispersion)$dispersion
        residual.scale <- as.vector(sqrt(dispersion))
        pred <- predict.lm(object, newdata, se.fit, scale = residual.scale,
            type = ifelse(type == "link", "response", type),
            terms = terms, na.action = na.action)
        fit <- pred$fit
        se.fit <- pred$se.fit
        switch(type, response = {
            se.fit <- se.fit * abs(family(object)$mu.eta(fit))
            fit <- family(object)$linkinv(fit)
        }, link = , terms = )
        if (missing(newdata) && !is.null(na.act)) {
            fit <- napredict(na.act, fit)
            se.fit <- napredict(na.act, se.fit)
        }
        pred <- list(fit = fit, se.fit = se.fit, residual.scale = residual.scale)
    }
    pred
}
<bytecode: 0x86468c0>
<environment: namespace:stats>


################################################################################ plot.lm
> getAnywhere(plot.lm)
A single object matching ‘plot.lm’ was found
It was found in the following places
  registered S3 method for plot from namespace stats
  namespace:stats
with value

function (x, which = c(1L:3L, 5L), caption = list("Residuals vs Fitted",
    "Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage",
    expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))),
    panel = if (add.smooth) panel.smooth else points, sub.caption = NULL,
    main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(),
    ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
    qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"),
    label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25)
{
    dropInf <- function(x, h) {
        if (any(isInf <- h >= 1)) {
            warning(gettextf("not plotting observations with leverage one:\n  %s",
                paste(which(isInf), collapse = ", ")), call. = FALSE,
                domain = NA)
            x[isInf] <- NaN
        }
        x
    }
    if (!inherits(x, "lm"))
        stop("use only with \"lm\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 6))
        stop("'which' must be in 1:6")
    isGlm <- inherits(x, "glm")
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    r <- residuals(x)
    yh <- predict(x)
    w <- weights(x)
    if (!is.null(w)) {
        wind <- w != 0
        r <- r[wind]
        yh <- yh[wind]
        w <- w[wind]
        labels.id <- labels.id[wind]
    }
    n <- length(r)
    if (any(show[2L:6L])) {
        s <- if (inherits(x, "rlm"))
            x$s
        else if (isGlm)
            sqrt(summary(x)$dispersion)
        else sqrt(deviance(x)/df.residual(x))
        hii <- (infl <- influence(x, do.coef = FALSE))$hat
        if (any(show[4L:6L])) {
            cook <- if (isGlm)
                cooks.distance(x, infl = infl)
            else cooks.distance(x, infl = infl, sd = s, res = r,
                hat = hii)
        }
    }
    if (any(show[2L:3L])) {
        ylab23 <- if (isGlm)
            "Std. deviance resid."
        else "Standardized residuals"
        r.w <- if (is.null(w))
            r
        else sqrt(w) * r
        rs <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
    }
    if (any(show[5L:6L])) {
        r.hat <- range(hii, na.rm = TRUE)
        isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 *
            mean(hii, na.rm = TRUE)
    }
    if (any(show[c(1L, 3L)]))
        l.fit <- if (isGlm)
            "Predicted values"
        else "Fitted values"
    if (is.null(id.n))
        id.n <- 0
    else {
        id.n <- as.integer(id.n)
        if (id.n < 0L || id.n > n)
            stop(gettextf("'id.n' must be in {1,..,%d}", n),
                domain = NA)
    }
    if (id.n > 0L) {
        if (is.null(labels.id))
            labels.id <- paste(1L:n)
        iid <- 1L:id.n
        show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
        if (any(show[2L:3L]))
            show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
        text.id <- function(x, y, ind, adj.x = TRUE) {
            labpos <- if (adj.x)
                label.pos[1 + as.numeric(x > mean(range(x)))]
            else 3
            text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
                pos = labpos, offset = 0.25)
        }
    }
    getCaption <- function(k) if (length(caption) < k)
        NA_character_
    else as.graphicsAnnot(caption[[k]])
    if (is.null(sub.caption)) {
        cal <- x$call
        if (!is.na(m.f <- match("formula", names(cal)))) {
            cal <- cal[c(1, m.f)]
            names(cal)[2L] <- ""
        }
        cc <- deparse(cal, 80)
        nc <- nchar(cc[1L], "c")
        abbr <- length(cc) > 1 || nc > 75
        sub.caption <- if (abbr)
            paste(substr(cc[1L], 1L, min(75L, nc)), "...")
        else cc[1L]
    }
    one.fig <- prod(par("mfcol")) == 1
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    if (show[1L]) {
        ylim <- range(r, na.rm = TRUE)
        if (id.n > 0)
            ylim <- extendrange(r = ylim, f = 0.08)
        dev.hold()
        plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main,
            ylim = ylim, type = "n", ...)
        panel(yh, r, ...)
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(getCaption(1), 3, 0.25, cex = cex.caption)
        if (id.n > 0) {
            y.id <- r[show.r]
            y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
            text.id(yh[show.r], y.id, show.r)
        }
        abline(h = 0, lty = 3, col = "gray")
        dev.flush()
    }
    if (show[2L]) {
        ylim <- range(rs, na.rm = TRUE)
        ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
        dev.hold()
        qq <- qqnorm(rs, main = main, ylab = ylab23, ylim = ylim,
            ...)
        if (qqline)
            qqline(rs, lty = 3, col = "gray50")
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(getCaption(2), 3, 0.25, cex = cex.caption)
        if (id.n > 0)
            text.id(qq$x[show.rs], qq$y[show.rs], show.rs)
        dev.flush()
    }
    if (show[3L]) {
        sqrtabsr <- sqrt(abs(rs))
        ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
        yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
        yhn0 <- if (is.null(w))
            yh
        else yh[w != 0]
        dev.hold()
        plot(yhn0, sqrtabsr, xlab = l.fit, ylab = yl, main = main,
            ylim = ylim, type = "n", ...)
        panel(yhn0, sqrtabsr, ...)
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(getCaption(3), 3, 0.25, cex = cex.caption)
        if (id.n > 0)
            text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs)
        dev.flush()
    }
    if (show[4L]) {
        if (id.n > 0) {
            show.r <- order(-cook)[iid]
            ymx <- cook[show.r[1L]] * 1.075
        }
        else ymx <- max(cook, na.rm = TRUE)
        dev.hold()
        plot(cook, type = "h", ylim = c(0, ymx), main = main,
            xlab = "Obs. number", ylab = "Cook's distance", ...)
        if (one.fig)
            title(sub = sub.caption, ...)
        mtext(getCaption(4), 3, 0.25, cex = cex.caption)
        if (id.n > 0)
            text.id(show.r, cook[show.r], show.r, adj.x = FALSE)
        dev.flush()
    }
    if (show[5L]) {
        ylab5 <- if (isGlm)
            "Std. Pearson resid."
        else "Standardized residuals"
        r.w <- residuals(x, "pearson")
        if (!is.null(w))
            r.w <- r.w[wind]
        rsp <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
        ylim <- range(rsp, na.rm = TRUE)
        if (id.n > 0) {
            ylim <- extendrange(r = ylim, f = 0.08)
            show.rsp <- order(-cook)[iid]
        }
        do.plot <- TRUE
        if (isConst.hat) {
            if (missing(caption))
                caption[[5L]] <- "Constant Leverage:\n Residuals vs Factor Levels"
            aterms <- attributes(terms(x))
            dcl <- aterms$dataClasses[-aterms$response]
            facvars <- names(dcl)[dcl %in% c("factor", "ordered")]
            mf <- model.frame(x)[facvars]
            if (ncol(mf) > 0) {
                dm <- data.matrix(mf)
                nf <- length(nlev <- unlist(unname(lapply(x$xlevels,
                  length))))
                ff <- if (nf == 1)
                  1
                else rev(cumprod(c(1, nlev[nf:2])))
                facval <- (dm - 1) %*% ff
                xx <- facval
                dev.hold()
                plot(facval, rsp, xlim = c(-1/2, sum((nlev -
                  1) * ff) + 1/2), ylim = ylim, xaxt = "n", main = main,
                  xlab = "Factor Level Combinations", ylab = ylab5,
                  type = "n", ...)
                axis(1, at = ff[1L] * (1L:nlev[1L] - 1/2) - 1/2,
                  labels = x$xlevels[[1L]])
                mtext(paste(facvars[1L], ":"), side = 1, line = 0.25,
                  adj = -0.05)
                abline(v = ff[1L] * (0:nlev[1L]) - 1/2, col = "gray",
                  lty = "F4")
                panel(facval, rsp, ...)
                abline(h = 0, lty = 3, col = "gray")
                dev.flush()
            }
            else {
                message(gettextf("hat values (leverages) are all = %s\n and there are no factor predictors; no plot no. 5",
                  format(mean(r.hat))), domain = NA)
                frame()
                do.plot <- FALSE
            }
        }
        else {
            xx <- hii
            xx[xx >= 1] <- NA
            dev.hold()
            plot(xx, rsp, xlim = c(0, max(xx, na.rm = TRUE)),
                ylim = ylim, main = main, xlab = "Leverage",
                ylab = ylab5, type = "n", ...)
            panel(xx, rsp, ...)
            abline(h = 0, v = 0, lty = 3, col = "gray")
            if (one.fig)
                title(sub = sub.caption, ...)
            if (length(cook.levels)) {
                p <- x$rank
                usr <- par("usr")
                hh <- seq.int(min(r.hat[1L], r.hat[2L]/100),
                  usr[2L], length.out = 101)
                for (crit in cook.levels) {
                  cl.h <- sqrt(crit * p * (1 - hh)/hh)
                  lines(hh, cl.h, lty = 2, col = 2)
                  lines(hh, -cl.h, lty = 2, col = 2)
                }
                legend("bottomleft", legend = "Cook's distance",
                  lty = 2, col = 2, bty = "n")
                xmax <- min(0.99, usr[2L])
                ymult <- sqrt(p * (1 - xmax)/xmax)
                aty <- sqrt(cook.levels) * ymult
                axis(4, at = c(-rev(aty), aty), labels = paste(c(rev(cook.levels),
                  cook.levels)), mgp = c(0.25, 0.25, 0), las = 2,
                  tck = 0, cex.axis = cex.id, col.axis = 2)
            }
            dev.flush()
        }
        if (do.plot) {
            mtext(getCaption(5), 3, 0.25, cex = cex.caption)
            if (id.n > 0) {
                y.id <- rsp[show.rsp]
                y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
                text.id(xx[show.rsp], y.id, show.rsp)
            }
        }
    }
    if (show[6L]) {
        g <- dropInf(hii/(1 - hii), hii)
        ymx <- max(cook, na.rm = TRUE) * 1.025
        dev.hold()
        plot(g, cook, xlim = c(0, max(g, na.rm = TRUE)), ylim = c(0,
            ymx), main = main, ylab = "Cook's distance", xlab = expression("Leverage  " *
            h[ii]), xaxt = "n", type = "n", ...)
        panel(g, cook, ...)
        athat <- pretty(hii)
        axis(1, at = athat/(1 - athat), labels = paste(athat))
        if (one.fig)
            title(sub = sub.caption, ...)
        p <- x$rank
        bval <- pretty(sqrt(p * cook/g), 5)
        usr <- par("usr")
        xmax <- usr[2L]
        ymax <- usr[4L]
        for (i in seq_along(bval)) {
            bi2 <- bval[i]^2
            if (p * ymax > bi2 * xmax) {
                xi <- xmax + strwidth(" ")/3
                yi <- bi2 * xi/p
                abline(0, bi2, lty = 2)
                text(xi, yi, paste(bval[i]), adj = 0, xpd = TRUE)
            }
            else {
                yi <- ymax - 1.5 * strheight(" ")
                xi <- p * yi/bi2
                lines(c(0, xi), c(0, yi), lty = 2)
                text(xi, ymax - 0.8 * strheight(" "), paste(bval[i]),
                  adj = 0.5, xpd = TRUE)
            }
        }
        mtext(getCaption(6), 3, 0.25, cex = cex.caption)
        if (id.n > 0) {
            show.r <- order(-cook)[iid]
            text.id(g[show.r], cook[show.r], show.r)
        }
        dev.flush()
    }
    if (!one.fig && par("oma")[3L] >= 1)
        mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
    invisible()
}
<bytecode: 0x86a0e48>
<environment: namespace:stats>



####################################################################### plot.lm
  A single object matching ‘plot.lm’ was found
It was found in the following places
registered S3 method for plot from namespace stats
namespace:stats
with value

function (x, which = c(1L:3L, 5L), caption = list("Residuals vs Fitted",
                                                  "Normal Q-Q", "Scale-Location", "Cook's distance", "Residuals vs Leverage",
                                                  expression("Cook's dist vs Leverage  " * h[ii]/(1 - h[ii]))),
          panel = if (add.smooth) panel.smooth else points, sub.caption = NULL,
          main = "", ask = prod(par("mfcol")) < length(which) && dev.interactive(),
          ..., id.n = 3, labels.id = names(residuals(x)), cex.id = 0.75,
          qqline = TRUE, cook.levels = c(0.5, 1), add.smooth = getOption("add.smooth"),
          label.pos = c(4, 2), cex.caption = 1, cex.oma.main = 1.25)
{
  dropInf <- function(x, h) {
    if (any(isInf <- h >= 1)) {
      warning(gettextf("not plotting observations with leverage one:\n  %s",
                       paste(which(isInf), collapse = ", ")), call. = FALSE,
              domain = NA)
      x[isInf] <- NaN
    }
    x
  }
  if (!inherits(x, "lm"))
    stop("use only with \"lm\" objects")
  if (!is.numeric(which) || any(which < 1) || any(which > 6))
    stop("'which' must be in 1:6")
  isGlm <- inherits(x, "glm")
  show <- rep(FALSE, 6)
  show[which] <- TRUE
  r <- residuals(x)
  yh <- predict(x)
  w <- weights(x)
  if (!is.null(w)) {
    wind <- w != 0
    r <- r[wind]
    yh <- yh[wind]
    w <- w[wind]
    labels.id <- labels.id[wind]
  }
  n <- length(r)
  if (any(show[2L:6L])) {
    s <- if (inherits(x, "rlm"))
      x$s
    else if (isGlm)
      sqrt(summary(x)$dispersion)
    else sqrt(deviance(x)/df.residual(x))
    hii <- (infl <- influence(x, do.coef = FALSE))$hat
    if (any(show[4L:6L])) {
      cook <- if (isGlm)
        cooks.distance(x, infl = infl)
      else cooks.distance(x, infl = infl, sd = s, res = r,
                          hat = hii)
    }
  }
  if (any(show[2L:3L])) {
    ylab23 <- if (isGlm)
      "Std. deviance resid."
    else "Standardized residuals"
    r.w <- if (is.null(w))
      r
    else sqrt(w) * r
    rs <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
  }
  if (any(show[5L:6L])) {
    r.hat <- range(hii, na.rm = TRUE)
    isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 *
      mean(hii, na.rm = TRUE)
  }
  if (any(show[c(1L, 3L)]))
    l.fit <- if (isGlm)
      "Predicted values"
  else "Fitted values"
  if (is.null(id.n))
    id.n <- 0
  else {
    id.n <- as.integer(id.n)
    if (id.n < 0L || id.n > n)
      stop(gettextf("'id.n' must be in {1,..,%d}", n),
           domain = NA)
  }
  if (id.n > 0L) {
    if (is.null(labels.id))
      labels.id <- paste(1L:n)
    iid <- 1L:id.n
    show.r <- sort.list(abs(r), decreasing = TRUE)[iid]
    if (any(show[2L:3L]))
      show.rs <- sort.list(abs(rs), decreasing = TRUE)[iid]
    text.id <- function(x, y, ind, adj.x = TRUE) {
      labpos <- if (adj.x)
        label.pos[1 + as.numeric(x > mean(range(x)))]
      else 3
      text(x, y, labels.id[ind], cex = cex.id, xpd = TRUE,
           pos = labpos, offset = 0.25)
    }
  }
  getCaption <- function(k) if (length(caption) < k)
    NA_character_
  else as.graphicsAnnot(caption[[k]])
  if (is.null(sub.caption)) {
    cal <- x$call
    if (!is.na(m.f <- match("formula", names(cal)))) {
      cal <- cal[c(1, m.f)]
      names(cal)[2L] <- ""
    }
    cc <- deparse(cal, 80)
    nc <- nchar(cc[1L], "c")
    abbr <- length(cc) > 1 || nc > 75
    sub.caption <- if (abbr)
      paste(substr(cc[1L], 1L, min(75L, nc)), "...")
    else cc[1L]
  }
  one.fig <- prod(par("mfcol")) == 1
  if (ask) {
    oask <- devAskNewPage(TRUE)
    on.exit(devAskNewPage(oask))
  }
  if (show[1L]) {
    ylim <- range(r, na.rm = TRUE)
    if (id.n > 0)
      ylim <- extendrange(r = ylim, f = 0.08)
    dev.hold()
    plot(yh, r, xlab = l.fit, ylab = "Residuals", main = main,
         ylim = ylim, type = "n", ...)
    panel(yh, r, ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(1), 3, 0.25, cex = cex.caption)
    if (id.n > 0) {
      y.id <- r[show.r]
      y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
      text.id(yh[show.r], y.id, show.r)
    }
    abline(h = 0, lty = 3, col = "gray")
    dev.flush()
  }
  if (show[2L]) {
    ylim <- range(rs, na.rm = TRUE)
    ylim[2L] <- ylim[2L] + diff(ylim) * 0.075
    dev.hold()
    qq <- qqnorm(rs, main = main, ylab = ylab23, ylim = ylim,
                 ...)
    if (qqline)
      qqline(rs, lty = 3, col = "gray50")
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(2), 3, 0.25, cex = cex.caption)
    if (id.n > 0)
      text.id(qq$x[show.rs], qq$y[show.rs], show.rs)
    dev.flush()
  }
  if (show[3L]) {
    sqrtabsr <- sqrt(abs(rs))
    ylim <- c(0, max(sqrtabsr, na.rm = TRUE))
    yl <- as.expression(substitute(sqrt(abs(YL)), list(YL = as.name(ylab23))))
    yhn0 <- if (is.null(w))
      yh
    else yh[w != 0]
    dev.hold()
    plot(yhn0, sqrtabsr, xlab = l.fit, ylab = yl, main = main,
         ylim = ylim, type = "n", ...)
    panel(yhn0, sqrtabsr, ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(3), 3, 0.25, cex = cex.caption)
    if (id.n > 0)
      text.id(yhn0[show.rs], sqrtabsr[show.rs], show.rs)
    dev.flush()
  }
  if (show[4L]) {
    if (id.n > 0) {
      show.r <- order(-cook)[iid]
      ymx <- cook[show.r[1L]] * 1.075
    }
    else ymx <- max(cook, na.rm = TRUE)
    dev.hold()
    plot(cook, type = "h", ylim = c(0, ymx), main = main,
         xlab = "Obs. number", ylab = "Cook's distance", ...)
    if (one.fig)
      title(sub = sub.caption, ...)
    mtext(getCaption(4), 3, 0.25, cex = cex.caption)
    if (id.n > 0)
      text.id(show.r, cook[show.r], show.r, adj.x = FALSE)
    dev.flush()
  }
  if (show[5L]) {
    ylab5 <- if (isGlm)
      "Std. Pearson resid."
    else "Standardized residuals"
    r.w <- residuals(x, "pearson")
    if (!is.null(w))
      r.w <- r.w[wind]
    rsp <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
    ylim <- range(rsp, na.rm = TRUE)
    if (id.n > 0) {
      ylim <- extendrange(r = ylim, f = 0.08)
      show.rsp <- order(-cook)[iid]
    }
    do.plot <- TRUE
    if (isConst.hat) {
      if (missing(caption))
        caption[[5L]] <- "Constant Leverage:\n Residuals vs Factor Levels"
      aterms <- attributes(terms(x))
      dcl <- aterms$dataClasses[-aterms$response]
      facvars <- names(dcl)[dcl %in% c("factor", "ordered")]
      mf <- model.frame(x)[facvars]
      if (ncol(mf) > 0) {
        dm <- data.matrix(mf)
        nf <- length(nlev <- unlist(unname(lapply(x$xlevels,
                                                  length))))
        ff <- if (nf == 1)
          1
        else rev(cumprod(c(1, nlev[nf:2])))
        facval <- (dm - 1) %*% ff
        xx <- facval
        dev.hold()
        plot(facval, rsp, xlim = c(-1/2, sum((nlev -
                                                1) * ff) + 1/2), ylim = ylim, xaxt = "n", main = main,
             xlab = "Factor Level Combinations", ylab = ylab5,
             type = "n", ...)
        axis(1, at = ff[1L] * (1L:nlev[1L] - 1/2) - 1/2,
             labels = x$xlevels[[1L]])
        mtext(paste(facvars[1L], ":"), side = 1, line = 0.25,
              adj = -0.05)
        abline(v = ff[1L] * (0:nlev[1L]) - 1/2, col = "gray",
               lty = "F4")
        panel(facval, rsp, ...)
        abline(h = 0, lty = 3, col = "gray")
        dev.flush()
      }
      else {
        message(gettextf("hat values (leverages) are all = %s\n and there are no factor predictors; no plot no. 5",
                         format(mean(r.hat))), domain = NA)
        frame()
        do.plot <- FALSE
      }
    }
    else {
      xx <- hii
      xx[xx >= 1] <- NA
      dev.hold()
      plot(xx, rsp, xlim = c(0, max(xx, na.rm = TRUE)),
           ylim = ylim, main = main, xlab = "Leverage",
           ylab = ylab5, type = "n", ...)
      panel(xx, rsp, ...)
      abline(h = 0, v = 0, lty = 3, col = "gray")
      if (one.fig)
        title(sub = sub.caption, ...)
      if (length(cook.levels)) {
        p <- x$rank
        usr <- par("usr")
        hh <- seq.int(min(r.hat[1L], r.hat[2L]/100),
                      usr[2L], length.out = 101)
        for (crit in cook.levels) {
          cl.h <- sqrt(crit * p * (1 - hh)/hh)
          lines(hh, cl.h, lty = 2, col = 2)
          lines(hh, -cl.h, lty = 2, col = 2)
        }
        legend("bottomleft", legend = "Cook's distance",
               lty = 2, col = 2, bty = "n")
        xmax <- min(0.99, usr[2L])
        ymult <- sqrt(p * (1 - xmax)/xmax)
        aty <- sqrt(cook.levels) * ymult
        axis(4, at = c(-rev(aty), aty), labels = paste(c(rev(cook.levels),
                                                         cook.levels)), mgp = c(0.25, 0.25, 0), las = 2,
             tck = 0, cex.axis = cex.id, col.axis = 2)
      }
      dev.flush()
    }
    if (do.plot) {
      mtext(getCaption(5), 3, 0.25, cex = cex.caption)
      if (id.n > 0) {
        y.id <- rsp[show.rsp]
        y.id[y.id < 0] <- y.id[y.id < 0] - strheight(" ")/3
        text.id(xx[show.rsp], y.id, show.rsp)
      }
    }
  }
  if (show[6L]) {
    g <- dropInf(hii/(1 - hii), hii)
    ymx <- max(cook, na.rm = TRUE) * 1.025
    dev.hold()
    plot(g, cook, xlim = c(0, max(g, na.rm = TRUE)), ylim = c(0,
                                                              ymx), main = main, ylab = "Cook's distance", xlab = expression("Leverage  " *
                                                                                                                               h[ii]), xaxt = "n", type = "n", ...)
    panel(g, cook, ...)
    athat <- pretty(hii)
    axis(1, at = athat/(1 - athat), labels = paste(athat))
    if (one.fig)
      title(sub = sub.caption, ...)
    p <- x$rank
    bval <- pretty(sqrt(p * cook/g), 5)
    usr <- par("usr")
    xmax <- usr[2L]
    ymax <- usr[4L]
    for (i in seq_along(bval)) {
      bi2 <- bval[i]^2
      if (p * ymax > bi2 * xmax) {
        xi <- xmax + strwidth(" ")/3
        yi <- bi2 * xi/p
        abline(0, bi2, lty = 2)
        text(xi, yi, paste(bval[i]), adj = 0, xpd = TRUE)
      }
      else {
        yi <- ymax - 1.5 * strheight(" ")
        xi <- p * yi/bi2
        lines(c(0, xi), c(0, yi), lty = 2)
        text(xi, ymax - 0.8 * strheight(" "), paste(bval[i]),
             adj = 0.5, xpd = TRUE)
      }
    }
    mtext(getCaption(6), 3, 0.25, cex = cex.caption)
    if (id.n > 0) {
      show.r <- order(-cook)[iid]
      text.id(g[show.r], cook[show.r], show.r)
    }
    dev.flush()
  }
  if (!one.fig && par("oma")[3L] >= 1)
    mtext(sub.caption, outer = TRUE, cex = cex.oma.main)
  invisible()
}

###################################################################################################
object=X<-glm(D~G*E,data=dat, family=binomial)
mat=list(D=dat$D[1:20], G=dat$G[1:20,], E=as.matrix(dat$E[1:20]))
predict(object, newdata, type="response", se=T)$se.fit
predict(object, newdata, type="link", se=T)$se.fit
predict(object, newdata=mat, type="link", interval="confidence", se.fit=T)
sqrt(diag(X %*% vcov(object) %*% t(X)))
plogis((diag(X %*% vcov(object) %*% t(X))))
###################################################################################################

thedata=read.csv("~/Research/Ray689/Salinity_Data_Modified_Leverage_Points.csv")
salinity  = thedata$salinity
lagged    = thedata$lagged.salinity
trend     = thedata$trend
discharge = thedata$discharge
n         = length(salinity)
ord       = order(discharge)
x         = discharge[ord]
y         = salinity[ord]
discharge = x
salinity  = y


fitlm = lm(salinity ~ bs(discharge,5))
t = qt(0.975,df=n-6)

predlm = predict(fitlm, se=T)
t(pred_vals <- model.matrix(fitlm) %*% coef(fitlm))
predlm$fit

t(resids <- salinity-pred_vals)
resid(fitlm)

# (sigmahat2 = sum((resids-mean(resids))^2)/(n-1))
var(resids)
(sigmahat2 = sum((resids-mean(resids))^2)/(n-6))

# XtXi = solve(t(model.matrix(fitlm)) %*% model.matrix(fitlm))

(se_fit = sqrt(diag(model.matrix(fitlm) %*% vcov(fitlm) %*% t(model.matrix(fitlm)))))
predlm$se.fit

(confint <- cbind(fit=pred_vals, lwr=pred_vals-t*se_fit, upr=pred_vals+t*se_fit))
(cf = predict(fitlm, interval="confidence"))
round((((cf[,3]-cf[,2])/2)/t)^2-predlm$se.fit^2, 14)


se_pred = sqrt((se_fit^2) + sigmahat2)

(predint <- cbind(fit=pred_vals, lwr=pred_vals-t*se_pred, upr=pred_vals+t*se_pred))
(pf = predict(fitlm, interval="prediction"))
(((pf[,3]-pf[,2])/2)/t)^2-se_fit^2
###########################################################################################################
###########################################################################################################
ppp = function(object, newdata, se.fit=FALSE,
                         interval=c("none", "confidence"), level=0.95,
                         type=c("response", "link"), na.action=na.pass, ...) {
  if(level>=1 & level<100) {level=level/100}
  stopifnot(0<level, level<1)
  if(pmatch(interval[[1]], "prediction", nomatch=0)) {stop("prediction intervals for binary response are not meaningful")}
  interval <- match.arg(interval)
  type <- match.arg(type)
  print(interval)
  print(type)
  if(se.fit && type=="response") {
    warning("standard errors can only be reported in the link scale")
    se.fit = FALSE
  }
  print(se.fit)
}

}}}}
}#}}
##############################################################################################
f=function(x1, x2, x3, x4, ...) {
cl <- match.call()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula", "data", "subset", "weights", "na.action", "offset"), names(mf), 0L)
return(list(cl, mf, m))
mf <- mf[c(1L, m)]
mf$drop.unused.levels <- TRUE
}

f=function(x1, x2, x3, x4, ...) {
  args(match.fun(match.call()[[1]]))
}

f=function(x1, x2=2, x3, x4){
  cl = match.call()
  nf = names(formals())
  m = match(nf, names(cl), -1L)
  print(list(name.cl=names(cl), cl=cl, nf=nf, m=m))
  undefined = nf[which(m==-1)]
  # for(i in seq_along(undefined)) assign(undefined[i], 77.7)
  lapply(nf[which(m==-1)], function(x) assign(x, NULL, pos=parent.frame(n=2)))
  cat(x1, x2, x3, x4)
}
f()

f=function(x1, x2=9, x3, x4){
  nf = names(formals())
  # for(i in seq_along(nf)) print(get0(nf[i], ifnotfound="missed"))
  # return(lapply(nf, function(x) get0(x, ifnotfound="missed")))
  # cat(x1, x2, x3, x4)
  print(enn(as.expression("x1")))
  print(enn("x2"))
  print(enn("x3"))
  print(enn("x4"))
  sapply(nf, enn)
}
f()
f(x3=4)
f(2)
f(1,2)
f(x2="TEE")
f(1,2,3,4)

cat("\nlala\n")
  lapply(nf[which(match(nf, names(match.call()), -1L) == -1)], function(x) assign(x, NULL, pos=parent.frame(n=2)))
  cat(x1, x2, x3, x4)
  is.null(x2)
}
f()
f(x3=4)
f(2)
f(1,2)
f(x2="TEE")
f(1,2,3,4)

# for(i in seq_along(nf)) print(get0(nf[i], ifnotfound="missed"))
# return(lapply(nf, function(x) get0(x, ifnotfound="missed")))
# get0("x3", ifnotfound="missed", envir=environment())

f=function(x1, x2=9, x3, x4){
  nf = names(formals())
  exists("x4", where=parent.frame(), inherits=FALSE)
}
f()
f(x3=4)
f(2)
f(1,2)
f(x2="TEE")
f(1,2,3,4)
1
f=function(x1, x2=9, x3, x4) {
  print(missing(x2))
  print(enn(x2))
  x2
}

f()
f(1)
f(x2=NULL)


f=function(x1, x2=9, x3, x4){
  nf = names(formals())
  lapply(nf[sapply(mget(nf), is.symbol)], function(x) assign(x, NULL, pos=parent.frame(n=2)))
  mget(nf)
  # cat(x1, x2, x3, x4)
}
f()
f("first")
f(1,2)
f(x3=4)
f(2)
f(1,2)
f(x2="TEE")
f(1,2,3,4)

a=function(...) f(...)
###############################################################################################

set.seed(2018)
dat = simulateCC(ncase=500, ncontrol=500, beta0=-4.165,
                 betaG_SNP=c(log(1.2), log(1.2), 0, log(1.2), 0),
                 betaE_bin=log(1.5),
                 betaGE_SNP_bin=c(log(1.3), 0, 0, log(1.3), 0),
                 MAF=c(0.1, 0.3, 0.3, 0.3, 0.1),
                 SNP_cor=0.7, E_bin_freq=0.5)

# SPMLE with known population disease rate of 0.03 and asymptotic SE estimates
spmleCombo(D=D, G=G, E=E, pi1=0.03, data=dat, nboot=0)

# Simulation with a single SNP and a single binary environmental variable.
# True population disease rate in this simulation is 0.03.
# This simulation scenario was used in the Supplementary Material of Stalder et. al. (2017)
# to compare performance against the less flexible method of Chatterjee and Carroll (2005),
# which is available as the function as snp.logistic in the Bioconductor package CGEN.
dat2 = simulateCC(ncase=500, ncontrol=500, beta0=-3.77,
                  betaG_SNP=log(1.2), betaE_bin=c(log(1.5), 0),
                  betaGE_SNP_bin=c(log(1.3), 0), MAF=0.1,
                  E_bin_freq=0.5)

# SPMLE using the rare disease assumption, 50 bootstraps, 2 cores
spmle(D=D, G=G, E=E, pi1=0, data=dat2)
spmleCombo(D=D, G=G, E=E, pi1=0, data=dat2, nboot=50, ncores=2)
spmleCombo(D=D, G=G, E=E, pi1=0.03, data=dat2, nboot=100, ncores=2)

########################################################################################
pack <- "caseControlGE"
path <- find.package(pack)
system(paste(shQuote(file.path(R.home("bin"), "R")), "CMD", "Rd2pdf", shQuote(path)))


