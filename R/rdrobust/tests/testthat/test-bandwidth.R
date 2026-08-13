## rdbwselect() structural contracts. Numerical values are pinned separately in
## test-numerical-contracts.R.

test_that("selected bandwidths are positive and finite", {
  fx <- make_rd_fixture(n = 1000)
  for (bw in c("mserd", "msetwo", "msesum", "msecomb1", "msecomb2",
               "cerrd", "certwo", "cersum", "cercomb1", "cercomb2")) {
    b <- rdbwselect(y = fx$y, x = fx$x, bwselect = bw)
    expect_true(all(is.finite(b$bws)), label = bw)
    expect_true(all(b$bws > 0), label = bw)
  }
})

test_that("all = TRUE returns every selector, matching the individual calls", {
  fx <- make_rd_fixture(n = 1000)
  bw_all <- rdbwselect(y = fx$y, x = fx$x, all = TRUE)
  expect_equal(nrow(bw_all$bws), 10L)

  for (bw in c("mserd", "msetwo", "cerrd")) {
    one <- rdbwselect(y = fx$y, x = fx$x, bwselect = bw)
    expect_equal(unname(bw_all$bws[rownames(bw_all$bws) == bw, , drop = TRUE]),
                 unname(one$bws[1, , drop = TRUE]),
                 tolerance = 1e-9, label = bw)
  }
})

test_that("msetwo allows the two sides to differ while mserd does not", {
  fx <- make_rd_fixture(n = 1000)
  rd <- rdbwselect(y = fx$y, x = fx$x, bwselect = "mserd")$bws
  two <- rdbwselect(y = fx$y, x = fx$x, bwselect = "msetwo")$bws
  expect_equal(unname(rd[1, 1]), unname(rd[1, 2]))
  expect_false(isTRUE(all.equal(unname(two[1, 1]), unname(two[1, 2]))))
})

test_that("CER bandwidths are smaller than the corresponding MSE bandwidths", {
  fx <- make_rd_fixture(n = 1000)
  mse <- rdbwselect(y = fx$y, x = fx$x, bwselect = "mserd")$bws[1, 1]
  cer <- rdbwselect(y = fx$y, x = fx$x, bwselect = "cerrd")$bws[1, 1]
  expect_lt(cer, mse)
})

test_that("rdrobust's internal selection agrees with a standalone rdbwselect", {
  fx <- make_rd_fixture(n = 1000)
  for (bw in c("mserd", "msetwo", "cerrd")) {
    est <- rdrobust(y = fx$y, x = fx$x, bwselect = bw)
    sel <- rdbwselect(y = fx$y, x = fx$x, bwselect = bw)
    expect_equal(unname(est$bws[1, ]), unname(sel$bws[1, 1:2]),
                 tolerance = 1e-10, label = bw)
    expect_equal(unname(est$bws[2, ]), unname(sel$bws[1, 3:4]),
                 tolerance = 1e-10, label = bw)
  }
})

test_that("a manual h is honoured and reported unchanged", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, h = 0.25)
  expect_equal(unname(est$bws[1, ]), c(0.25, 0.25))
  expect_equal(est$bwselect, "Manual")
})

test_that("a manual h with b is honoured on both rows", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, h = 0.25, b = 0.4)
  expect_equal(unname(est$bws[1, ]), c(0.25, 0.25))
  expect_equal(unname(est$bws[2, ]), c(0.4, 0.4))
})

test_that("rho sets b as h/rho", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, rho = 0.5)
  expect_equal(unname(est$bws[2, 1]), 0.6, tolerance = 1e-12)
})

test_that("a two-sided h vector is applied side by side", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, h = c(0.2, 0.35))
  expect_equal(unname(est$bws[1, ]), c(0.2, 0.35))
  expect_equal(unname(est$N_h[1]),
               sum(fx$x < fx$c & abs(fx$x - fx$c) <= 0.2))
  expect_equal(unname(est$N_h[2]),
               sum(fx$x >= fx$c & abs(fx$x - fx$c) <= 0.35))
})

test_that("bwrestrict caps the bandwidth at the support of x", {
  fx <- make_rd_fixture(n = 1000)
  b <- rdbwselect(y = fx$y, x = fx$x, bwselect = "mserd", bwrestrict = TRUE)
  bw_max <- max(abs(fx$c - min(fx$x)), abs(fx$c - max(fx$x)))
  expect_lte(max(b$bws), bw_max)
})
