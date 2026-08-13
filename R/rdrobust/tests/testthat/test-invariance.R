## Metamorphic contracts: transformations of the input that must leave the
## output unchanged (or change it in an exactly known way). These catch whole
## classes of indexing, sorting and scaling defects that fixed-value baselines
## cannot.

test_that("results do not depend on the row order of the data", {
  fx <- make_rd_fixture(n = 1000)
  set.seed(101)
  o <- sample(fx$n)

  expect_equal(rd_key(rdrobust(y = fx$y[o], x = fx$x[o])),
               rd_key(rdrobust(y = fx$y, x = fx$x)), tolerance = 1e-10)

  expect_equal(rd_key(rdrobust(y = fx$y[o], x = fx$x[o], covs = fx$Z[o, ])),
               rd_key(rdrobust(y = fx$y, x = fx$x, covs = fx$Z)),
               tolerance = 1e-10)

  expect_equal(rd_key(rdrobust(y = fx$y[o], x = fx$x[o],
                               cluster = fx$cluster[o], vce = "cr1")),
               rd_key(rdrobust(y = fx$y, x = fx$x,
                               cluster = fx$cluster, vce = "cr1")),
               tolerance = 1e-10)

  expect_equal(rd_key(rdrobust(y = fx$y[o], x = fx$x[o],
                               weights = fx$weights[o])),
               rd_key(rdrobust(y = fx$y, x = fx$x, weights = fx$weights)),
               tolerance = 1e-10)

  expect_equal(rd_key(rdrobust(y = fx$y[o], x = fx$x[o], vce = "hc3")),
               rd_key(rdrobust(y = fx$y, x = fx$x, vce = "hc3")),
               tolerance = 1e-10)
})

test_that("subset= is equivalent to subsetting the inputs by hand", {
  fx <- make_rd_fixture(n = 1000)
  s <- fx$x > -0.6

  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x, subset = s)),
               rd_key(rdrobust(y = fx$y[s], x = fx$x[s])), tolerance = 1e-10)

  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x, covs = fx$Z, subset = s)),
               rd_key(rdrobust(y = fx$y[s], x = fx$x[s], covs = fx$Z[s, ])),
               tolerance = 1e-10)

  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x, cluster = fx$cluster,
                               vce = "cr1", subset = s)),
               rd_key(rdrobust(y = fx$y[s], x = fx$x[s],
                               cluster = fx$cluster[s], vce = "cr1")),
               tolerance = 1e-10)
})

test_that("the estimate is equivariant in the scale and location of y", {
  fx <- make_rd_fixture(n = 1000)
  base <- rd_key(rdrobust(y = fx$y, x = fx$x))

  expect_equal(rd_key(rdrobust(y = 3 * fx$y, x = fx$x)), 3 * base,
               tolerance = 1e-9)
  expect_equal(rd_key(rdrobust(y = fx$y + 10, x = fx$x)), base,
               tolerance = 1e-9)
})

test_that("the estimate is invariant to the scale and location of x", {
  fx <- make_rd_fixture(n = 1000)
  base <- rd_key(rdrobust(y = fx$y, x = fx$x))

  expect_equal(rd_key(rdrobust(y = fx$y, x = 5 * fx$x, c = 0)), base,
               tolerance = 1e-9)
  expect_equal(rd_key(rdrobust(y = fx$y, x = fx$x + 2, c = 2)), base,
               tolerance = 1e-9)
})

test_that("extreme magnitudes of x and y do not destabilise the estimate", {
  ## Guards the stdvars standardization default.
  fx <- make_rd_fixture(n = 1000)
  base <- rd_key(rdrobust(y = fx$y, x = fx$x))

  for (s in c(1e-6, 1e6)) {
    expect_equal(rd_key(rdrobust(y = s * fx$y, x = fx$x)), s * base,
                 tolerance = 1e-7, label = sprintf("y * %g", s))
    expect_equal(rd_key(rdrobust(y = fx$y, x = s * fx$x, c = 0)), base,
                 tolerance = 1e-7, label = sprintf("x * %g", s))
  }
})
