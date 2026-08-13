## Input validation. Every check that currently works is pinned here so that a
## refactor cannot silently remove it; the checks that are still missing live in
## test-known-bugs.R.

test_that("y and x must have the same length", {
  fx <- make_rd_fixture(n = 200)
  expect_error(rdrobust(y = fx$y[-1], x = fx$x))
})

test_that("the cutoff must be inside the support of x", {
  fx <- make_rd_fixture(n = 500)
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, c = 99)))
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, c = -99)))
})

test_that("the cutoff must be a single finite number", {
  fx <- make_rd_fixture(n = 500)
  expect_error(rdrobust(y = fx$y, x = fx$x, c = c(0, 1)))
  expect_error(rdrobust(y = fx$y, x = fx$x, c = NA))
  expect_error(rdrobust(y = fx$y, x = fx$x, c = Inf))
})

test_that("p, q and deriv are validated", {
  fx <- make_rd_fixture(n = 500)
  expect_error(rdrobust(y = fx$y, x = fx$x, p = -1))
  expect_error(rdrobust(y = fx$y, x = fx$x, deriv = -1))
  ## q must not be below p.
  expect_error(rdrobust(y = fx$y, x = fx$x, p = 2, q = 1))
  ## deriv must not exceed p.
  expect_error(rdrobust(y = fx$y, x = fx$x, p = 1, deriv = 2))
})

test_that("an unknown kernel is rejected", {
  fx <- make_rd_fixture(n = 500)
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, kernel = "bogus")))
})

test_that("an unknown bwselect is rejected", {
  fx <- make_rd_fixture(n = 500)
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, bwselect = "bogus")))
  expect_error(suppressWarnings(rdbwselect(y = fx$y, x = fx$x, bwselect = "bogus")))
})

test_that("an unknown vce is rejected", {
  fx <- make_rd_fixture(n = 500)
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, vce = "bogus")))
})

test_that("the confidence level must be inside (0, 100)", {
  fx <- make_rd_fixture(n = 500)
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, level = 0)))
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, level = 100)))
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, level = -5)))
})

test_that("a non-positive bandwidth is rejected", {
  fx <- make_rd_fixture(n = 500)
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, h = 0)))
  expect_error(suppressWarnings(rdrobust(y = fx$y, x = fx$x, h = -0.3)))
})

test_that("kernel and bwselect are case-insensitive in rdrobust", {
  fx <- make_rd_fixture(n = 500)
  a <- rdrobust(y = fx$y, x = fx$x, h = 0.3, kernel = "TRI")
  b <- rdrobust(y = fx$y, x = fx$x, h = 0.3, kernel = "tri")
  expect_equal(rd_key(a), rd_key(b), tolerance = 1e-12)
})

test_that("missing values in y or x are dropped consistently", {
  fx <- make_rd_fixture(n = 800)
  y <- fx$y; x <- fx$x
  y[1:10] <- NA
  x[11:20] <- NA
  keep <- !is.na(y) & !is.na(x)

  est <- rdrobust(y = y, x = x, h = 0.3)
  ref <- rdrobust(y = y[keep], x = x[keep], h = 0.3)
  expect_equal(rd_key(est), rd_key(ref), tolerance = 1e-12)
  expect_equal(sum(unname(est$N)), sum(keep))
})

test_that("missing values in covs drop the whole row", {
  fx <- make_rd_fixture(n = 800)
  Z <- fx$Z
  Z[1:15, 1] <- NA
  keep <- complete.cases(Z)

  est <- rdrobust(y = fx$y, x = fx$x, covs = Z, h = 0.3)
  ref <- rdrobust(y = fx$y[keep], x = fx$x[keep], covs = Z[keep, ], h = 0.3)
  expect_equal(rd_key(est), rd_key(ref), tolerance = 1e-12)
})

test_that("a covs matrix with the wrong number of rows is rejected", {
  fx <- make_rd_fixture(n = 500)
  expect_error(rdrobust(y = fx$y, x = fx$x, covs = fx$Z[-1, ]))
})

test_that("a cluster variable of the wrong length is rejected", {
  fx <- make_rd_fixture(n = 500)
  expect_error(rdrobust(y = fx$y, x = fx$x, cluster = fx$cluster[-1]))
})
