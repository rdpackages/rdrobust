## rdrobust() point estimates must equal an independently computed weighted
## least squares fit on each side of the cutoff. This is the primary
## correctness contract: it does not rely on any rdrobust internal.

test_that("sharp RD point estimates match the WLS oracle across kernel, p, deriv", {
  fx <- make_rd_fixture()
  h <- 0.3

  for (kernel in c("uni", "tri", "epa")) {
    for (p in 1:3) {
      for (deriv in 0:min(2, p - 1)) {
        est <- rdrobust(y = fx$y, x = fx$x, h = h, p = p, q = p + 1,
                        deriv = deriv, kernel = kernel)
        expect_equal(
          unname(est$coef[1]),
          rd_oracle(fx$y, fx$x, fx$c, h, p, deriv, kernel),
          tolerance = 1e-9,
          label = sprintf("kernel=%s p=%d deriv=%d", kernel, p, deriv)
        )
      }
    }
  }
})

test_that("sharp RD with covariates matches the stacked common-gamma oracle", {
  fx <- make_rd_fixture()
  h <- 0.3

  ## deriv >= 2 is excluded here on purpose: it is broken upstream.
  ## See test-known-bugs.R (XL-1).
  for (p in 1:3) {
    for (deriv in 0:min(1, p - 1)) {
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

test_that("fuzzy RD point estimate is the ratio of the two sharp jumps", {
  fz <- make_fuzzy_fixture()
  h <- 0.3

  for (p in 1:2) {
    est <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, h = h, p = p, q = p + 1)
    expect_equal(
      unname(est$coef[1]),
      rd_oracle_fuzzy(fz$y, fz$x, fz$d, 0, h, p, 0, "tri"),
      tolerance = 1e-9,
      label = sprintf("fuzzy p=%d", p)
    )
  }
})

test_that("a non-zero cutoff shifts nothing but the cutoff", {
  fx <- make_rd_fixture()
  shifted <- rdrobust(y = fx$y, x = fx$x + 5, c = 5, h = 0.3)
  base <- rdrobust(y = fx$y, x = fx$x, c = 0, h = 0.3)
  expect_equal(unname(shifted$coef[1]), unname(base$coef[1]), tolerance = 1e-10)
  expect_equal(unname(shifted$se[3]), unname(base$se[3]), tolerance = 1e-10)
})

test_that("scalepar rescales the point estimate linearly", {
  fx <- make_rd_fixture()
  base <- rdrobust(y = fx$y, x = fx$x, h = 0.3)
  for (s in c(-1, 0.5, 3)) {
    est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, scalepar = s)
    expect_equal(unname(est$coef[1]), s * unname(base$coef[1]),
                 tolerance = 1e-10, label = sprintf("scalepar=%g", s))
  }
})

test_that("left and right effective sample sizes count the kernel support", {
  fx <- make_rd_fixture()
  h <- 0.3
  est <- rdrobust(y = fx$y, x = fx$x, h = h)
  expect_equal(unname(est$N_h[1]), sum(fx$x < fx$c & abs(fx$x - fx$c) <= h))
  expect_equal(unname(est$N_h[2]), sum(fx$x >= fx$c & abs(fx$x - fx$c) <= h))
})
