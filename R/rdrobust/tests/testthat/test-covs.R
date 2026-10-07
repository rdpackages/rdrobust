## Covariate-adjustment contracts: interfaces, collinearity handling, and the
## effect on the estimate.

test_that("covs accepts a matrix, a data frame, a formula and a character vector", {
  fx <- make_rd_fixture(n = 1000)
  df <- data.frame(y = fx$y, x = fx$x, z1 = fx$Z[, 1], z2 = fx$Z[, 2])

  m <- rdrobust(y = fx$y, x = fx$x, h = 0.3, covs = fx$Z)
  d <- rdrobust(y = fx$y, x = fx$x, h = 0.3, covs = as.data.frame(fx$Z))
  f <- rdrobust(y = df$y, x = df$x, h = 0.3, covs = ~ z1 + z2, data = df)
  s <- rdrobust(y = df$y, x = df$x, h = 0.3, covs = c("z1", "z2"), data = df)

  expect_equal(rd_key(d), rd_key(m), tolerance = 1e-10)
  expect_equal(rd_key(f), rd_key(m), tolerance = 1e-10)
  expect_equal(rd_key(s), rd_key(m), tolerance = 1e-10)
})

test_that("a formula expands factors and transformations", {
  fx <- make_rd_fixture(n = 1000)
  g <- factor(rep_len(c("a", "b", "c"), fx$n))
  df <- data.frame(y = fx$y, x = fx$x, z1 = fx$Z[, 1], g = g)

  f <- rdrobust(y = df$y, x = df$x, h = 0.3, covs = ~ z1 + g, data = df)
  m <- rdrobust(y = fx$y, x = fx$x, h = 0.3,
                covs = cbind(fx$Z[, 1], model.matrix(~ g)[, -1]))
  expect_equal(rd_key(f), rd_key(m), tolerance = 1e-10)
})

test_that("adding a relevant covariate changes the estimate", {
  fx <- make_rd_fixture(n = 1000)
  plain <- rdrobust(y = fx$y, x = fx$x, h = 0.3)
  adj <- rdrobust(y = fx$y, x = fx$x, h = 0.3, covs = fx$Z)
  expect_false(isTRUE(all.equal(unname(plain$coef[1]), unname(adj$coef[1]))))
})

test_that("a covariate orthogonal to everything barely moves the estimate", {
  fx <- make_rd_fixture(n = 1000)
  plain <- rdrobust(y = fx$y, x = fx$x, h = 0.3)
  noise <- rdrobust(y = fx$y, x = fx$x, h = 0.3,
                    covs = matrix(fx$Z[, 2], ncol = 1))
  expect_equal(unname(noise$coef[1]), unname(plain$coef[1]), tolerance = 0.1)
})

test_that("collinear covariates are dropped with a warning by default", {
  fx <- make_rd_fixture(n = 1000)
  a <- fx$Z[, 1]
  expect_warning(
    dup <- rdrobust(y = fx$y, x = fx$x, h = 0.3, covs = cbind(a, a)),
    "ollinear"
  )
  single <- rdrobust(y = fx$y, x = fx$x, h = 0.3, covs = matrix(a, ncol = 1))
  expect_equal(unname(dup$coef[1]), unname(single$coef[1]), tolerance = 1e-8)
})

test_that("a constant covariate column is absorbed without error", {
  fx <- make_rd_fixture(n = 1000)
  est <- suppressWarnings(
    rdrobust(y = fx$y, x = fx$x, h = 0.3,
             covs = cbind(fx$Z[, 1], rep(1, fx$n)))
  )
  expect_true(is.finite(est$coef[1]))
})

test_that("the covariate-adjusted model is labelled as such", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, covs = fx$Z)
  expect_true(grepl("ovariate", est$rdmodel))
})

test_that("covariates are passed through to bandwidth selection", {
  fx <- make_rd_fixture(n = 1000)
  h_c <- rdbwselect(y = fx$y, x = fx$x, covs = fx$Z)$bws[1, 1]
  h_0 <- rdbwselect(y = fx$y, x = fx$x)$bws[1, 1]
  expect_false(isTRUE(all.equal(h_c, h_0)))
})
