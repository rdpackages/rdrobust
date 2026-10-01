## Open defects, encoded as executable specifications.
##
## Each test below asserts the CORRECT behaviour and is skipped with the audit
## id of the defect it guards. When a defect is fixed, delete the skip() line:
## the test then either passes (fix confirmed) or fails (fix incomplete).
##
## Reference: AUDIT_FINDINGS_rdrobust_2026-07-20.md and
## AUDIT_STATUS_2026-08-12.md in the project working directory.

test_that("XL-1: factorial(deriv) is applied in the sharp+covariates branch", {

  fx <- make_rd_fixture()
  h <- 0.3
  for (p in 3:4) {
    for (deriv in 2:min(3, p - 1)) {
      est <- rdrobust(y = fx$y, x = fx$x, covs = fx$Z, h = h, p = p,
                      q = p + 1, deriv = deriv, kernel = "tri")
      expect_equal(
        unname(est$coef[1]),
        rd_oracle_covs(fx$y, fx$x, fx$c, fx$Z, h, p, deriv, "tri"),
        tolerance = 1e-9,
        label = sprintf("covs p=%d deriv=%d", p, deriv)
      )
    }
  }
})

test_that("R-3: sharpbw with a non-zero cutoff does not crash", {

  fz <- make_fuzzy_fixture()
  expect_no_error(
    rdrobust(y = fz$y, x = fz$x + 5, c = 5, fuzzy = fz$d, sharpbw = TRUE)
  )
})

test_that("R-5: vce is lowercased before the mass-point pre-sort decision", {

  fx <- make_rd_fixture()
  lower <- rdrobust(y = fx$y, x = fx$x, vce = "nn", masspoints = "off")
  upper <- rdrobust(y = fx$y, x = fx$x, vce = "NN", masspoints = "off")
  expect_equal(rd_key(upper), rd_key(lower), tolerance = 1e-10)
})

test_that("R/PY-2a: rdbwselect reports the effective sample size, not the full sample", {

  fx <- make_rd_fixture(c = 1000)
  fx$x <- 1000 + 250 * (fx$x - 1000)
  bw <- rdbwselect(y = fx$y, x = fx$x, c = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, c = 1000)
  expect_equal(unname(bw$N_h), unname(est$N_h))
})

test_that("R/PY-2b: rdbwselect returns the cutoff on the original scale", {

  fx <- make_rd_fixture()
  fx$x <- 1000 + 250 * fx$x
  bw <- rdbwselect(y = fx$y, x = fx$x, c = 1000)
  expect_equal(bw$c, 1000)
})

test_that("R-9: rdbwselect accepts the same option casing as rdrobust", {

  fx <- make_rd_fixture()
  expect_no_error(rdbwselect(y = fx$y, x = fx$x, kernel = "TRI"))
  expect_equal(
    rdbwselect(y = fx$y, x = fx$x, kernel = "TRI")$bws[1, 1],
    rdbwselect(y = fx$y, x = fx$x, kernel = "tri")$bws[1, 1]
  )
})

test_that("NEW-1: a degenerate bandwidth is reported, never returned as NaN", {

  ## Two values of x per side cannot identify the pilot polynomials. The
  ## selector must either return a finite bandwidth or say so; never NaN, and
  ## never a crash from deep inside the fit.
  mp <- make_masspoint_fixture(n = 2000, K = 5, seed = 2)
  bw <- tryCatch(suppressWarnings(rdbwselect(y = mp$y, x = mp$x)),
                 error = function(e) conditionMessage(e))
  if (is.character(bw)) expect_match(bw, "Not enough variability")
  else expect_true(is.finite(bw$bws[1, 1]))
})

test_that("NEW-2: covs_drop=FALSE with collinear covariates errors informatively", {

  fx <- make_rd_fixture()
  a <- rnorm(fx$n)
  expect_error(
    rdrobust(y = fx$y, x = fx$x, covs = cbind(a, a), covs_drop = FALSE, h = 0.3),
    "collinear|multicollinear|rank"
  )
})

test_that("NEW-5: a bandwidth that empties one side errors informatively", {

  mp <- make_masspoint_fixture(n = 1000, K = 5, seed = 3)
  expect_error(rdrobust(y = mp$y, x = mp$x, h = 0.3),
               "observations|side|bandwidth")
})

test_that("R-8: vcov.rdrobust documents its diagonal-only contract", {
  ## RESOLVED by documenting, which was the audit's stated alternative to
  ## filling: the design matrices needed to form the cross-covariances are not
  ## retained on the fitted object. The zeros are now labelled as placeholders
  ## so they cannot be mistaken for estimated independence.
  fx <- make_rd_fixture()
  est <- rdrobust(y = fx$y, x = fx$x, h = 0.3)
  V <- vcov(est)
  expect_equal(diag(V), unname(as.vector(est$se))^2, tolerance = 1e-12,
               ignore_attr = TRUE)
  expect_true(all(V[upper.tri(V)] == 0))
  expect_match(attr(V, "offdiag"), "not estimated")
})

test_that("R-10: nnmatch is validated", {

  ## rdrobust's established pattern: the specific complaint is a warning, and
  ## the abort itself is the generic "invalid input" stop().
  fx <- make_rd_fixture()
  for (bad in list(0, -3, 2.5, c(1, 2))) {
    expect_warning(
      try(rdrobust(y = fx$y, x = fx$x, h = 0.3, nnmatch = bad), silent = TRUE),
      "nnmatch", label = toString(bad)
    )
    expect_error(
      suppressWarnings(rdrobust(y = fx$y, x = fx$x, h = 0.3, nnmatch = bad)),
      "invalid input", label = toString(bad)
    )
  }
})

test_that("R-10: negative weights are rejected rather than silently dropped", {

  fx <- make_rd_fixture()
  w <- fx$weights
  w[1:50] <- -1
  expect_error(rdrobust(y = fx$y, x = fx$x, h = 0.3, weights = w), "weight")
})

test_that("R-12: rdplot validates binselect", {

  fx <- make_rd_fixture()
  expect_error(rdplot(y = fx$y, x = fx$x, binselect = "bogus", hide = TRUE),
               "binselect")
})
