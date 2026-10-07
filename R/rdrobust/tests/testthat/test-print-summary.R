## Print and summary methods, and the accessor S3 methods.

test_that("print and summary run for all three classes", {
  fx <- make_rd_fixture(n = 800)
  est <- rdrobust(y = fx$y, x = fx$x)
  bw <- rdbwselect(y = fx$y, x = fx$x)
  pl <- rdplot(y = fx$y, x = fx$x, hide = TRUE)

  expect_output(print(est))
  expect_output(print(summary(est)))
  expect_output(print(bw))
  expect_output(print(summary(bw)))
  expect_output(print(pl))
  expect_output(print(summary(pl)))
})

test_that("summary(all = TRUE) prints the three inference rows", {
  fx <- make_rd_fixture(n = 800)
  est <- rdrobust(y = fx$y, x = fx$x)

  compact <- capture.output(print(summary(est)))
  full <- capture.output(print(summary(est, all = TRUE)))
  detail <- capture.output(print(summary(est, detail = TRUE)))

  ## Compact: a single "RD Effect" row, no method breakdown.
  expect_true(any(grepl("RD Effect", compact)))
  expect_false(any(grepl("Conventional", compact)))

  ## detail = TRUE: Conventional + Robust, but not Bias-Corrected.
  expect_true(any(grepl("Conventional", detail)))
  expect_true(any(grepl("Robust", detail)))
  expect_false(any(grepl("Bias-Corrected|Bias-corrected", detail)))

  ## all = TRUE: all three rows.
  expect_true(any(grepl("Conventional", full)))
  expect_true(any(grepl("Bias-Corrected|Bias-corrected", full)))
  expect_true(any(grepl("Robust", full)))
  expect_gt(length(full), length(detail))
})

test_that("summary output does not depend on the estimate being reprinted", {
  fx <- make_rd_fixture(n = 800)
  est <- rdrobust(y = fx$y, x = fx$x)
  a <- capture.output(print(summary(est)))
  b <- capture.output(print(summary(est)))
  expect_identical(a, b)
})

test_that("coef and vcov accessors return the documented shapes", {
  fx <- make_rd_fixture(n = 800)
  est <- rdrobust(y = fx$y, x = fx$x)

  cf <- coef(est)
  V <- vcov(est)
  expect_true(is.numeric(cf))
  expect_true(is.matrix(V))
  expect_equal(nrow(V), ncol(V))
  expect_true(all(diag(V) > 0))
})

test_that("the fitted object exposes the documented slots", {
  fx <- make_rd_fixture(n = 800)
  est <- rdrobust(y = fx$y, x = fx$x)
  for (nm in c("coef", "se", "z", "pv", "ci", "bws", "N", "N_h", "c", "p", "q",
               "kernel", "bwselect", "vce", "rdmodel")) {
    expect_true(!is.null(est[[nm]]), label = nm)
  }
  expect_equal(dim(est$ci), c(3L, 2L))
  expect_equal(dim(est$bws), c(2L, 2L))
})

test_that("rdbwselect exposes the documented slots", {
  fx <- make_rd_fixture(n = 800)
  bw <- rdbwselect(y = fx$y, x = fx$x)
  for (nm in c("bws", "bwselect", "kernel", "p", "q", "c", "N", "vce")) {
    expect_true(!is.null(bw[[nm]]), label = nm)
  }
  expect_equal(ncol(bw$bws), 4L)
})
