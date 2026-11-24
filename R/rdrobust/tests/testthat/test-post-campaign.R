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

## PR #25 review (Matias, 2026-10-01).

test_that("Review 3: the few-cluster warning ignores zero-weight observations", {

  i <- seq_len(800)
  x <- (i - 400.5) / 400
  y <- 1 + x + 0.7 * (x >= 0) + 0.2 * sin(i * 0.113)
  cluster <- i %% 40
  weights <- as.numeric(cluster < 4)  # 4 contributing clusters per side
  expect_warning(
    rdrobust(y, x, h = 0.8, cluster = cluster, weights = weights, vce = "cr1"),
    "Only 4 (left) and 4 (right) clusters", fixed = TRUE
  )
})

test_that("Review 1: b = h/rho for non-round rho, scalar and asymmetric h", {

  fx <- make_rd_fixture()
  for (rho in c(0.123456, 0.00014, 0.00001, 0.5)) {
    est <- rdrobust(y = fx$y, x = fx$x, h = 0.5, rho = rho)
    expect_equal(unname(est$bws[2, ]), rep(0.5 / rho, 2), tolerance = 1e-12)
    est <- rdrobust(y = fx$y, x = fx$x, h = c(0.4, 0.6), rho = rho)
    expect_equal(unname(est$bws[2, ]), c(0.4, 0.6) / rho, tolerance = 1e-12)
  }
})

test_that("D3: p+1 mass-point clusters warn about a degenerate variance", {

  ## Two mass points per side within h (cluster = x), p = 1: the fit is
  ## saturated at the cluster level and the conventional SE is near zero.
  x <- rep(c(-5:-1, 0:4), each = 30)
  set.seed(3)
  y <- 1 + 0.5 * x + (x >= 0.5) + rnorm(length(x))
  w <- character()
  est <- withCallingHandlers(
    rdrobust(y, x, c = 0.5, h = 2, b = 5, p = 1, cluster = x, vce = "cr1"),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("Only 2 (left) and 2 (right) clusters", w, fixed = TRUE) &
                  grepl("estimated variance may be degenerate", w, fixed = TRUE)))
  expect_lt(est$se[1], 1e-8)
})

test_that("D3: p+1 general clusters can have positive variance", {
  x <- seq(-1, 1, length.out = 800)
  i <- 0:799
  cluster <- i %% 2 + 2 * (x >= 0)
  y <- 1 + x + 0.7 * (x >= 0) + 0.3 * (i %% 2) + 0.2 * sin(i * 0.113)
  w <- character()
  est <- withCallingHandlers(
    rdrobust(y, x, h = 0.8, b = 0.9, cluster = cluster, vce = "cr1"),
    warning = function(cnd) { w <<- c(w, conditionMessage(cnd)); invokeRestart("muffleWarning") }
  )
  expect_true(any(grepl("inference may be unreliable", w, fixed = TRUE)))
  expect_false(any(grepl("not identified", w, fixed = TRUE)))

  ## Independent local-linear CR1 calculation using the estimation sample.
  variances <- vapply(list(x < 0, x >= 0), function(side) {
    keep <- side & abs(x) < 0.8
    X <- cbind(1, x[keep])
    weights <- pmax(1 - abs(x[keep]) / 0.8, 0)
    invG <- solve(crossprod(X, weights * X))
    residual <- y[keep] - X %*% (invG %*% crossprod(X, weights * y[keep]))
    scores <- vapply(unique(cluster[keep]), function(g) {
      selected <- cluster[keep] == g
      as.vector(crossprod(X[selected, ], weights[selected] * residual[selected]))
    }, numeric(2))
    V <- invG %*% tcrossprod(scores) %*% invG * ((sum(keep) - 1) / (sum(keep) - 2)) * 2
    V[1, 1]
  }, numeric(1))
  expect_true(all(variances > 0))
  expect_equal(unname(est$se[1]), sqrt(sum(variances)), tolerance = 1e-12)
})

test_that("conventional HC1 and CR1 SEs use the estimation window", {
  i <- seq_len(1000)
  x <- (i - 500.5) / 100
  y <- 1 + 0.2 * x + 1.42 * (x >= 0) + sin(i / 3) + cos(i / 8)
  cluster <- (i - 1L) %/% 10L

  ## Reference from separate fits inside h using sandwich::vcovHC and
  ## sandwich::vcovCL, both with type = "HC1".
  expected <- c(hc1 = 0.259027485692751, cr1 = 0.711581929424620)
  for (vce in names(expected)) {
    ids <- if (vce == "cr1") cluster else NULL
    for (b in c(0.8, 1.2, 3.8)) {
      est <- rdrobust(y, x, h = 1.2, b = b, kernel = "uniform",
                      cluster = ids, vce = vce, masspoints = "off")
      expect_equal(est$se[1], unname(expected[vce]), tolerance = 1e-10)
    }

    ## The fuzzy first stage uses the same adjustment, also with asymmetric
    ## bandwidths and zero user weights.
    fuzzy <- as.numeric(i %% 11L < ifelse(x >= 0, 8L, 4L))
    weights <- as.numeric(i %% 7L != 0L)
    h <- c(1.2, 1.6)
    fit <- function(b) rdrobust(y, x, fuzzy = fuzzy, weights = weights,
                               h = h, b = b, kernel = "uniform", cluster = ids,
                               vce = vce, masspoints = "off")
    narrow <- fit(h)
    wide <- fit(c(3.2, 3.8))
    expect_equal(wide$se[1], narrow$se[1], tolerance = 1e-10)
    expect_equal(wide$se_T[1], narrow$se_T[1], tolerance = 1e-10)
  }
})

test_that("a wider bias window cannot supply conventional degrees of freedom", {
  x <- c(-10:-1, 1:10)
  y <- 1 + x + (x >= 0) + sin(x)
  hc1 <- rdrobust(y, x, h = 2.5, b = 7.5, kernel = "uniform", vce = "hc1")
  expect_true(is.nan(hc1$se[1]))
  expect_true(is.finite(hc1$se[3]))

  ids <- sign(x) * ceiling(abs(x) / 3)
  expect_warning(
    cr1 <- rdrobust(y, x, h = 3.5, b = 7.5, kernel = "uniform",
                    cluster = ids, vce = "cr1"),
    "estimated variance may be degenerate"
  )
  expect_true(is.nan(cr1$se[1]))
  expect_true(is.finite(cr1$se[3]))
})

test_that("Review 2: all = TRUE stops when any selector is undefined; mserd alone works", {

  set.seed(1)
  x <- sample(-7:6, 400, TRUE)
  y <- 0.3 * x + (x >= 0) + rnorm(400)
  expect_error(suppressWarnings(rdbwselect(y, x, p = 2, all = TRUE)), "Not enough variability")
  expect_error(suppressWarnings(rdbwselect(y, x, p = 2, bwselect = "msecomb1")), "Not enough variability")
  bw <- suppressWarnings(rdbwselect(y, x, p = 2))$bws
  expect_true(all(is.finite(bw) & bw > 0))
})

test_that("Review: a cutoff outside the support is reported as such", {

  fx <- make_rd_fixture()
  expect_error(rdbwselect(fx$y, fx$x, c = 2), "within the range of x")
  expect_warning(try(rdrobust(fx$y, fx$x, c = 2), silent = TRUE), "within the range of x")
})

test_that("Review: rdbwselect validates nnmatch and weights like rdrobust", {

  fx <- make_rd_fixture()
  for (bad in list(0, 2.5)) {
    expect_error(suppressWarnings(rdbwselect(fx$y, fx$x, nnmatch = bad)), "invalid input",
                 label = toString(bad))
  }
  w <- fx$weights; w[1:5] <- -1
  expect_error(rdbwselect(fx$y, fx$x, weights = w), "non-negative")
  expect_no_error(suppressWarnings(rdbwselect(fx$y, fx$x, masspoints = "")))
})
