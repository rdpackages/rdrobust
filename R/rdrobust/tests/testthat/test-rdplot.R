## rdplot() structural contracts. hide = TRUE keeps the test suite headless.

test_that("rdplot returns the documented structure", {
  fx <- make_rd_fixture(n = 800)
  pl <- rdplot(y = fx$y, x = fx$x, hide = TRUE)
  expect_s3_class(pl, "rdplot")
  expect_true(!is.null(pl$vars_bins))
  expect_true(!is.null(pl$vars_poly))
  expect_true(!is.null(pl$J))
})

test_that("nbins is honoured on both sides", {
  fx <- make_rd_fixture(n = 800)
  pl <- rdplot(y = fx$y, x = fx$x, nbins = c(12, 15), hide = TRUE)
  expect_equal(unname(pl$J), c(12, 15))
})

test_that("every bin mean lies inside the range of its own observations", {
  fx <- make_rd_fixture(n = 800)
  pl <- rdplot(y = fx$y, x = fx$x, nbins = c(10, 10), hide = TRUE)
  vb <- pl$vars_bins
  expect_true(all(vb$rdplot_mean_y >= vb$rdplot_min_y - 1e-10, na.rm = TRUE))
  expect_true(all(vb$rdplot_mean_y <= vb$rdplot_max_y + 1e-10, na.rm = TRUE))
  expect_true(all(vb$rdplot_mean_x >= vb$rdplot_min_x - 1e-10, na.rm = TRUE))
  expect_true(all(vb$rdplot_mean_x <= vb$rdplot_max_x + 1e-10, na.rm = TRUE))
})

test_that("bin observation counts account for every observation", {
  ## rdplot bins the full support: unlike rdrobust, h controls the polynomial
  ## fit, not which observations are binned.
  fx <- make_rd_fixture(n = 800)
  pl <- rdplot(y = fx$y, x = fx$x, nbins = c(10, 10), hide = TRUE)
  expect_equal(sum(pl$vars_bins$rdplot_N, na.rm = TRUE), fx$n)

  n_left <- sum(fx$x < fx$c)
  left_bins <- pl$vars_bins$rdplot_mean_x < fx$c
  expect_equal(sum(pl$vars_bins$rdplot_N[left_bins], na.rm = TRUE), n_left)
})

test_that("the documented binselect methods all run", {
  fx <- make_rd_fixture(n = 800)
  for (bs in c("es", "espr", "esmv", "esmvpr",
               "qs", "qspr", "qsmv", "qsmvpr")) {
    pl <- rdplot(y = fx$y, x = fx$x, binselect = bs, hide = TRUE)
    expect_true(all(pl$J > 0), label = bs)
  }
})

test_that("quantile-spaced bins hold roughly equal counts", {
  fx <- make_rd_fixture(n = 800)
  pl <- rdplot(y = fx$y, x = fx$x, binselect = "qs", nbins = c(10, 10),
               hide = TRUE)
  n_per_bin <- pl$vars_bins$rdplot_N
  n_per_bin <- n_per_bin[!is.na(n_per_bin) & n_per_bin > 0]
  expect_lt(max(n_per_bin) / min(n_per_bin), 2)
})

test_that("the polynomial fit is evaluated on both sides of the cutoff", {
  fx <- make_rd_fixture(n = 800)
  pl <- rdplot(y = fx$y, x = fx$x, p = 4, hide = TRUE)
  vp <- pl$vars_poly
  expect_true(any(vp$rdplot_x < fx$c))
  expect_true(any(vp$rdplot_x >= fx$c))
  expect_true(all(is.finite(vp$rdplot_y)))
})

test_that("the confidence-interval columns appear when ci is requested", {
  fx <- make_rd_fixture(n = 800)
  pl <- rdplot(y = fx$y, x = fx$x, ci = 95, hide = TRUE)
  expect_true(!is.null(pl$vars_bins$rdplot_ci_l))
  expect_true(!is.null(pl$vars_bins$rdplot_ci_r))
  ok <- !is.na(pl$vars_bins$rdplot_ci_l)
  expect_true(all(pl$vars_bins$rdplot_ci_l[ok] <=
                  pl$vars_bins$rdplot_ci_r[ok]))
})

test_that("rdplot respects subset=", {
  fx <- make_rd_fixture(n = 800)
  s <- fx$x > -0.7
  a <- rdplot(y = fx$y, x = fx$x, subset = s, nbins = c(8, 8), hide = TRUE)
  b <- rdplot(y = fx$y[s], x = fx$x[s], nbins = c(8, 8), hide = TRUE)
  expect_equal(a$vars_bins$rdplot_mean_y, b$vars_bins$rdplot_mean_y,
               tolerance = 1e-10)
})

test_that("bin edges stay attached to their own bin when a side has holes", {
  ## The left support has a gap in [-0.55, -0.25], so several evenly-spaced
  ## left bins come out empty while the right side is fully occupied. Empty
  ## bins are dropped from vars_bins, and the surviving rows must keep the
  ## edges of the bins they were actually computed from.
  set.seed(20260813)
  x <- c(runif(200, -1, -0.55), runif(200, -0.25, 0), runif(400, 0, 1))
  y <- 3 + 2 * x + 4 * (x >= 0) + rnorm(length(x), sd = 0.3)

  pl <- rdplot(y = y, x = x, nbins = c(12, 12), binselect = "es", hide = TRUE)
  vb <- pl$vars_bins

  ## The gap must actually have emptied some left bins, or the test is vacuous.
  expect_lt(sum(vb$rdplot_mean_x < 0), pl$J[1])

  ## Each bin's mean falls inside that bin, and the reported edges bracket the
  ## midpoint that rdplot computed independently.
  expect_true(all(vb$rdplot_mean_x >= vb$rdplot_min_bin))
  expect_true(all(vb$rdplot_mean_x <= vb$rdplot_max_bin))
  expect_equal((vb$rdplot_min_bin + vb$rdplot_max_bin) / 2,
               vb$rdplot_mean_bin, tolerance = 1e-12)
})
