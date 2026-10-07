## Mass-point handling on discrete running variables.

test_that("a discrete running variable triggers the mass-point warning", {
  mp <- make_masspoint_fixture(n = 800, K = 30)
  expect_warning(rdrobust(y = mp$y, x = mp$x), "[Mm]ass point")
})

test_that("masspoints='off' suppresses the mass-point warning", {
  mp <- make_masspoint_fixture(n = 800, K = 30)
  expect_silent(rdrobust(y = mp$y, x = mp$x, masspoints = "off"))
})

test_that("a continuous running variable produces no mass-point warning", {
  fx <- make_rd_fixture(n = 800)
  expect_silent(rdrobust(y = fx$y, x = fx$x))
})

test_that("masspoints='adjust' changes the selected bandwidth on discrete x", {
  mp <- make_masspoint_fixture(n = 800, K = 30)
  adj <- suppressWarnings(rdbwselect(y = mp$y, x = mp$x,
                                     masspoints = "adjust"))$bws[1, 1]
  off <- rdbwselect(y = mp$y, x = mp$x, masspoints = "off")$bws[1, 1]
  expect_false(isTRUE(all.equal(adj, off)))
})

test_that("the number of unique mass points is reported per side", {
  mp <- make_masspoint_fixture(n = 800, K = 30)
  est <- suppressWarnings(rdrobust(y = mp$y, x = mp$x))
  expect_equal(unname(est$M[1]), length(unique(mp$x[mp$x < 0])))
  expect_equal(unname(est$M[2]), length(unique(mp$x[mp$x >= 0])))
})

test_that("mass points do not disturb the point estimate at a fixed bandwidth", {
  mp <- make_masspoint_fixture(n = 800, K = 30)
  est <- suppressWarnings(rdrobust(y = mp$y, x = mp$x, h = 0.5))
  expect_equal(unname(est$coef[1]),
               rd_oracle(mp$y, mp$x, 0, 0.5, 1, 0, "tri"),
               tolerance = 1e-9)
})
