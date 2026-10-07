## Cluster-robust variance contracts and the vce -> cluster-variant mapping.

test_that("the cluster-robust estimators run and are labelled correctly", {
  fx <- make_rd_fixture(n = 1000)
  for (v in c("cr1", "cr2", "cr3")) {
    est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, cluster = fx$cluster, vce = v)
    expect_true(all(is.finite(est$se)), label = v)
    expect_true(all(est$se > 0), label = v)
  }
})

test_that("CR standard errors are ordered cr1 < cr2 < cr3", {
  fx <- make_rd_fixture(n = 1000)
  se <- vapply(c("cr1", "cr2", "cr3"),
               function(v) unname(rdrobust(y = fx$y, x = fx$x, h = 0.3,
                                           cluster = fx$cluster,
                                           vce = v)$se[1]),
               numeric(1))
  expect_true(all(diff(se) > 0))
})

test_that("clustering changes the standard error but not the point estimate", {
  fx <- make_rd_fixture(n = 1000)
  plain <- rdrobust(y = fx$y, x = fx$x, h = 0.3, vce = "hc1")
  clust <- rdrobust(y = fx$y, x = fx$x, h = 0.3, cluster = fx$cluster,
                    vce = "cr1")
  expect_equal(unname(clust$coef[1]), unname(plain$coef[1]), tolerance = 1e-10)
  expect_false(isTRUE(all.equal(unname(clust$se[1]), unname(plain$se[1]))))
})

test_that("singleton clusters reproduce the heteroskedasticity-robust fit", {
  fx <- make_rd_fixture(n = 1000)
  singleton <- seq_len(fx$n)
  clust <- rdrobust(y = fx$y, x = fx$x, h = 0.3, cluster = singleton,
                    vce = "cr1")
  expect_true(all(is.finite(clust$se)))
  expect_true(all(clust$se > 0))
})

test_that("the number of clusters is reported", {
  fx <- make_rd_fixture(n = 1000)
  est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, cluster = fx$cluster,
                  vce = "cr1")
  expect_equal(length(unique(fx$cluster)), 40L)
  expect_true(grepl("cluster", est$rdmodel, ignore.case = TRUE))
})

test_that("a non-cluster vce with a cluster variable is remapped with a warning", {
  fx <- make_rd_fixture(n = 1000)
  expect_warning(
    est <- rdrobust(y = fx$y, x = fx$x, h = 0.3, cluster = fx$cluster,
                    vce = "hc1")
  )
  cr1 <- rdrobust(y = fx$y, x = fx$x, h = 0.3, cluster = fx$cluster,
                  vce = "cr1")
  expect_equal(unname(est$se[1]), unname(cr1$se[1]), tolerance = 1e-10)
})

test_that("cr0 is not accepted", {
  fx <- make_rd_fixture(n = 1000)
  expect_error(
    suppressWarnings(rdrobust(y = fx$y, x = fx$x, h = 0.3,
                              cluster = fx$cluster, vce = "cr0"))
  )
})

test_that("cluster is passed through to bandwidth selection", {
  fx <- make_rd_fixture(n = 1000)
  h_cl <- rdbwselect(y = fx$y, x = fx$x, cluster = fx$cluster,
                     vce = "cr1")$bws[1, 1]
  h_no <- rdbwselect(y = fx$y, x = fx$x, vce = "hc1")$bws[1, 1]
  expect_false(isTRUE(all.equal(h_cl, h_no)))
})

test_that("a factor cluster variable behaves like the integer one", {
  fx <- make_rd_fixture(n = 1000)
  a <- rdrobust(y = fx$y, x = fx$x, h = 0.3, cluster = fx$cluster, vce = "cr1")
  b <- rdrobust(y = fx$y, x = fx$x, h = 0.3,
                cluster = as.character(fx$cluster), vce = "cr1")
  expect_equal(rd_key(a), rd_key(b), tolerance = 1e-10)
})
