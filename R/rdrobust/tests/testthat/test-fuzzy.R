## Fuzzy RD contracts.

test_that("perfect compliance reproduces the sharp estimate", {
  fx <- make_rd_fixture(n = 1000)
  d <- as.numeric(fx$x >= fx$c)
  sharp <- rdrobust(y = fx$y, x = fx$x, h = 0.3)
  fuzzy <- suppressWarnings(
    rdrobust(y = fx$y, x = fx$x, fuzzy = d, h = 0.3)
  )
  expect_equal(unname(fuzzy$coef[1]), unname(sharp$coef[1]), tolerance = 1e-8)
})

test_that("the fuzzy estimate scales inversely with the first stage", {
  fz <- make_fuzzy_fixture(n = 1000)
  est <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, h = 0.3)
  first <- rdrobust(y = fz$d, x = fz$x, h = 0.3)
  reduced <- rdrobust(y = fz$y, x = fz$x, h = 0.3)
  expect_equal(unname(est$coef[1]),
               unname(reduced$coef[1]) / unname(first$coef[1]),
               tolerance = 1e-8)
})

test_that("the fuzzy model is labelled and runs with covariates and clusters", {
  fz <- make_fuzzy_fixture(n = 1000)
  est <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, h = 0.3)
  expect_true(grepl("Fuzzy", est$rdmodel))

  Z <- matrix(rnorm(fz$n), ncol = 1)
  with_covs <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, covs = Z, h = 0.3)
  expect_true(all(is.finite(with_covs$se)))

  cl <- rep_len(seq_len(30), fz$n)
  with_cl <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, cluster = cl,
                      vce = "cr1", h = 0.3)
  expect_true(all(is.finite(with_cl$se)))
})

test_that("fuzzy bandwidth selection runs for both sharpbw settings at c = 0", {
  fz <- make_fuzzy_fixture(n = 1000)
  a <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, sharpbw = FALSE)
  b <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, sharpbw = TRUE)
  expect_true(all(is.finite(a$bws)))
  expect_true(all(is.finite(b$bws)))
  expect_false(isTRUE(all.equal(unname(a$bws[1, 1]), unname(b$bws[1, 1]))))
})

test_that("fuzzy kink RD runs at deriv = 1", {
  fz <- make_fuzzy_fixture(n = 1000)
  est <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, deriv = 1, p = 2, h = 0.4)
  expect_true(is.finite(est$coef[1]))
  expect_true(grepl("Kink", est$rdmodel))
})
