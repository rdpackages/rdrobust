## Shared fixtures, kernel helper, and independent local-polynomial oracles.
##
## The oracles below recompute the RD estimand from first principles with plain
## weighted least squares. They are deliberately independent of the package
## internals: if rdrobust() and the oracle agree, the estimator is right for
## reasons that do not depend on rdrobust's own code being right.

## Kernel weights on the support |u| <= 1. Scaling is irrelevant for point
## estimates (WLS is invariant to a common positive scale factor), so the 1/h
## normalization used internally by rdrobust_kweight() is omitted here.
rd_w_fun <- function(u, kernel) {
  switch(kernel,
         "uni" = 0.5 * (abs(u) <= 1),
         "tri" = (1 - abs(u)) * (abs(u) <= 1),
         "epa" = 0.75 * (1 - u^2) * (abs(u) <= 1))
}

## Sharp RD fixture. Treatment effect 0.7 at the cutoff, one relevant and one
## irrelevant covariate, a cluster id, and positive weights.
make_rd_fixture <- function(n = 1500, seed = 7, c = 0) {
  set.seed(seed)
  x <- runif(n, -1, 1) + c
  Z <- cbind(z1 = rnorm(n), z2 = runif(n))
  y <- 1 + 2 * (x - c) + (x >= c) * 0.7 + 0.4 * Z[, 1] + rnorm(n)
  list(y = y, x = x, c = c, Z = Z, n = n,
       cluster = rep_len(seq_len(40), n),
       weights = runif(n, 0.5, 1.5))
}

## Fuzzy RD fixture with imperfect compliance on both sides.
make_fuzzy_fixture <- function(n = 1500, seed = 11) {
  set.seed(seed)
  x <- runif(n, -1, 1)
  d <- as.numeric(runif(n) < ifelse(x >= 0, 0.85, 0.20))
  y <- 1 + 2 * x + 0.7 * d + rnorm(n)
  list(y = y, x = x, d = d, n = n)
}

## Discrete running variable with a controllable number of mass points.
make_masspoint_fixture <- function(n = 800, K = 30, seed = 5) {
  set.seed(seed)
  x <- sample(seq(-1, 1, length.out = K), n, replace = TRUE)
  y <- 1 + 2 * x + (x >= 0) * 0.7 + rnorm(n)
  list(y = y, x = x, n = n, K = K)
}

## One-sided local polynomial fit; returns the deriv-th coefficient
## (not yet multiplied by factorial(deriv)).
rd_one_side <- function(y, x, c, h, p, deriv, kernel, side) {
  sub <- if (side == "r") x >= c else x < c
  xx <- x[sub] - c
  yy <- y[sub]
  w <- rd_w_fun(xx / h, kernel)
  keep <- w > 0
  X <- cbind(1, poly(xx[keep], degree = p, raw = TRUE))
  W <- w[keep]
  beta <- solve(crossprod(X * sqrt(W)), crossprod(X * W, yy[keep]))
  beta[deriv + 1]
}

## Oracle for the sharp RD estimand without covariates.
rd_oracle <- function(y, x, c, h, p, deriv, kernel) {
  factorial(deriv) *
    (rd_one_side(y, x, c, h, p, deriv, kernel, "r") -
     rd_one_side(y, x, c, h, p, deriv, kernel, "l"))
}

## Oracle for the sharp RD estimand with covariates. rdrobust fits separate
## polynomials on each side but a COMMON covariate coefficient vector, which is
## exactly the stacked design built here.
rd_oracle_covs <- function(y, x, c, Z, h, p, deriv, kernel) {
  xx <- x - c
  Tr <- as.numeric(x >= c)
  w <- rd_w_fun(xx / h, kernel)
  keep <- w > 0
  P <- poly(xx[keep], degree = p, raw = TRUE)
  Tk <- Tr[keep]
  Dl <- cbind(1, P) * (1 - Tk)
  Dr <- cbind(1, P) * Tk
  X <- cbind(Dl, Dr, Z[keep, , drop = FALSE])
  W <- w[keep]
  beta <- solve(crossprod(X * sqrt(W)), crossprod(X * W, y[keep]))
  factorial(deriv) * (beta[(p + 1) + deriv + 1] - beta[deriv + 1])
}

## Oracle for the fuzzy RD estimand: ratio of sharp reduced-form to sharp
## first-stage jumps at the cutoff.
rd_oracle_fuzzy <- function(y, x, d, c, h, p, deriv, kernel) {
  rd_oracle(y, x, c, h, p, deriv, kernel) /
    rd_oracle(d, x, c, h, p, deriv, kernel)
}

## Compact signature of an rdrobust fit, for invariance comparisons.
rd_key <- function(m) {
  c(tau = unname(m$coef[1]), tau_bc = unname(m$coef[2]),
    se = unname(m$se[1]), se_rb = unname(m$se[3]))
}
