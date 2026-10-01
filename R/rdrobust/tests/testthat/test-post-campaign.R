## Defects found after the 2026-08 fix campaign (D1-D3, R-17).

test_that("D2: a covariate that is near-collinear on one side does not make h depend on noise", {

  ## z2 equals z1 up to 1e-9 noise on the left only, so the pooled collinearity
  ## check keeps it while the left-side residualized Gram is numerically
  ## singular. With ginv.tol = 1e-20 that noise direction was inverted and the
  ## selected bandwidth moved with the noise draw (by about 6e-2).
  set.seed(1)
  n <- 2000
  x <- runif(n, -1, 1); z1 <- rnorm(n); zr <- rnorm(n)
  y <- 1 + x + 0.5 * (x >= 0) + 0.3 * z1 + rnorm(n)
  hs <- sapply(1:4, function(s) {
    set.seed(s)
    z2 <- ifelse(x < 0, z1 + 1e-9 * rnorm(n), zr)
    rdbwselect(y, x, covs = cbind(z1, z2))$bws[1, 1]
  })
  expect_lt(diff(range(hs)), 1e-6)
})

test_that("D3: few clusters within the bandwidth trigger a warning", {

  fx <- make_rd_fixture()
  few <- as.integer(cut(fx$x, breaks = seq(-1, 1, length.out = 9)))  # 4 per side
  expect_warning(
    rdrobust(y = fx$y, x = fx$x, h = 0.5, cluster = few, vce = "cr1"),
    "clusters"
  )
  ## The fixture's 40 clusters per side do not warn.
  expect_no_warning(
    rdrobust(y = fx$y, x = fx$x, h = 0.5, cluster = fx$cluster, vce = "cr1")
  )
})

test_that("D3: a cluster-robust standard error of exactly zero is flagged", {

  fx <- make_rd_fixture()
  two <- as.integer(cut(fx$x, breaks = c(-1, -0.5, 0, 0.5, 1), include.lowest = TRUE))
  w <- character(0)
  est <- withCallingHandlers(
    rdrobust(y = fx$y, x = fx$x, h = 1, p = 1, cluster = two, vce = "cr1"),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  if (any(est$se == 0)) expect_true(any(grepl("exactly 0", w)))
  expect_true(any(grepl("clusters", w)))
})

test_that("R-17: a fuzzy side with a single observation fails informatively", {

  fz <- make_fuzzy_fixture()
  keep <- fz$x >= 0 | seq_along(fz$x) == which(fz$x < 0)[1]
  err <- tryCatch(
    suppressWarnings(rdrobust(y = fz$y[keep], x = fz$x[keep], fuzzy = fz$d[keep])),
    error = function(e) conditionMessage(e)
  )
  expect_type(err, "character")
  expect_false(grepl("missing value where TRUE/FALSE needed", err))
})

test_that("Very discrete running variables give a bandwidth or an informative error", {

  ## With a few mass points the pilot polynomials are not identified. This used
  ## to end in NaN, in svd() errors, or in "missing value where TRUE/FALSE
  ## needed" from the nearest-neighbour loop.
  for (K in 3:9) for (seed in 1:2) for (p in 1:2) {
    mp <- make_masspoint_fixture(n = 2000, K = K, seed = seed)
    out <- tryCatch(suppressWarnings(rdbwselect(mp$y, mp$x, p = p)$bws[1, ]),
                    error = function(e) conditionMessage(e))
    lab <- sprintf("K=%d seed=%d p=%d", K, seed, p)
    if (is.character(out)) expect_match(out, "Not enough", label = lab)
    else expect_true(all(is.finite(out) & out > 0), label = lab)
  }
})
