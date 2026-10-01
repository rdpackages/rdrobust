## Weight handling contracts.

test_that("unit weights reproduce the unweighted fit exactly", {
  fx <- make_rd_fixture(n = 1000)
  w1 <- rep(1, fx$n)

  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x, weights = w1)),
               rd_key(rdrobust(y = fx$y, x = fx$x)), tolerance = 1e-12)

  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x, covs = fx$Z, weights = w1)),
               rd_key(rdrobust(y = fx$y, x = fx$x, covs = fx$Z)),
               tolerance = 1e-12)
})

test_that("a constant weight cancels out of the estimate", {
  fx <- make_rd_fixture(n = 1000)
  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x, weights = rep(2, fx$n))),
               rd_key(rdrobust(y = fx$y, x = fx$x)), tolerance = 1e-9)
})

test_that("weights reach the point estimate through the kernel weights", {
  fx <- make_rd_fixture(n = 1000)
  h <- 0.3
  est <- rdrobust(y = fx$y, x = fx$x, h = h, weights = fx$weights)

  ## Oracle: user weights multiply the kernel weights.
  side <- function(s) {
    sub <- if (s == "r") fx$x >= fx$c else fx$x < fx$c
    xx <- fx$x[sub] - fx$c
    w <- rd_w_fun(xx / h, "tri") * fx$weights[sub]
    keep <- w > 0
    X <- cbind(1, xx[keep])
    beta <- solve(crossprod(X * sqrt(w[keep])),
                  crossprod(X * w[keep], fx$y[sub][keep]))
    beta[1]
  }
  expect_equal(unname(est$coef[1]), side("r") - side("l"), tolerance = 1e-9)
})

test_that("weights change the selected bandwidth", {
  fx <- make_rd_fixture(n = 1000)
  hw <- rdbwselect(y = fx$y, x = fx$x, weights = fx$weights)$bws[1, 1]
  h0 <- rdbwselect(y = fx$y, x = fx$x)$bws[1, 1]
  expect_false(isTRUE(all.equal(hw, h0)))
})

test_that("zero weights drop observations from the fit", {
  fx <- make_rd_fixture(n = 1000)
  w <- rep(1, fx$n)
  drop <- fx$x > 0.1 & fx$x < 0.2
  w[drop] <- 0

  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x, h = 0.3, weights = w)),
               rd_key(rdrobust(y = fx$y[!drop], x = fx$x[!drop], h = 0.3)),
               tolerance = 1e-9)
})
