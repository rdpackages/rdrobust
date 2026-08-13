## Variance estimation contracts for the non-clustered estimators.

test_that("all documented vce options run and give positive standard errors", {
  fx <- make_rd_fixture(n = 1000)
  for (v in c("nn", "hc0", "hc1", "hc2", "hc3")) {
    est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, vce = v)
    expect_true(all(is.finite(est$se)), label = v)
    expect_true(all(est$se > 0), label = v)
  }
})

test_that("HC standard errors are ordered hc0 < hc1 < hc2 < hc3", {
  fx <- make_rd_fixture(n = 1000)
  se <- vapply(c("hc0", "hc1", "hc2", "hc3"),
               function(v) unname(rdrobust(y = fx$y, x = fx$x, h = 0.3,
                                           vce = v)$se[1]),
               numeric(1))
  expect_true(all(diff(se) > 0))
})

test_that("the conventional HC0 variance matches a direct sandwich computation", {
  fx <- make_rd_fixture(n = 1000)
  h <- 0.3
  est <- rdrobust(y = fx$y, x = fx$x, h = h, p = 1, vce = "hc0", kernel = "uni")

  ## Sum of the two one-sided sandwich variances for the intercept.
  side_var <- function(s) {
    sub <- if (s == "r") fx$x >= fx$c else fx$x < fx$c
    xx <- fx$x[sub] - fx$c
    yy <- fx$y[sub]
    w <- rd_w_fun(xx / h, "uni")
    keep <- w > 0
    X <- cbind(1, xx[keep]); W <- diag(w[keep]); yk <- yy[keep]
    XtWX_inv <- solve(t(X) %*% W %*% X)
    e <- as.numeric(yk - X %*% (XtWX_inv %*% t(X) %*% W %*% yk))
    meat <- t(X) %*% (W %*% diag(e^2) %*% W) %*% X
    (XtWX_inv %*% meat %*% XtWX_inv)[1, 1]
  }
  expect_equal(unname(est$se[1]), sqrt(side_var("l") + side_var("r")),
               tolerance = 1e-8)
})

test_that("the robust bias-corrected SE differs from the conventional one", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, b = 0.5)
  expect_false(isTRUE(all.equal(unname(est$se[1]), unname(est$se[3]))))
  expect_true(all(est$se > 0))
})

test_that("the three inference rows are mutually consistent", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x)

  ## Conventional and bias-corrected share a standard error; the estimates differ.
  expect_equal(unname(est$se[1]), unname(est$se[2]), tolerance = 1e-12)
  expect_equal(unname(est$coef[2]), unname(est$coef[3]), tolerance = 1e-12)

  ## z statistics and p-values are consistent with coef and se.
  expect_equal(unname(est$z), unname(est$coef / est$se), tolerance = 1e-10)
  expect_equal(unname(est$pv), unname(2 * pnorm(-abs(est$z))), tolerance = 1e-10)
})

test_that("confidence intervals match the reported level", {
  fx <- make_rd_fixture(n = 1000)
  for (lev in c(90, 95, 99)) {
    est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, level = lev)
    z <- qnorm(1 - (1 - lev / 100) / 2)
    expect_equal(unname(est$ci[3, 1]),
                 unname(est$coef[3] - z * est$se[3]), tolerance = 1e-9,
                 label = sprintf("level=%d lower", lev))
    expect_equal(unname(est$ci[3, 2]),
                 unname(est$coef[3] + z * est$se[3]), tolerance = 1e-9,
                 label = sprintf("level=%d upper", lev))
  }
})

test_that("nnmatch changes the nearest-neighbour variance", {
  fx <- make_rd_fixture(n = 1000)
  se3 <- unname(rdrobust(y = fx$y, x = fx$x, h = 0.3, vce = "nn",
                         nnmatch = 3)$se[1])
  se8 <- unname(rdrobust(y = fx$y, x = fx$x, h = 0.3, vce = "nn",
                         nnmatch = 8)$se[1])
  expect_false(isTRUE(all.equal(se3, se8)))
})
