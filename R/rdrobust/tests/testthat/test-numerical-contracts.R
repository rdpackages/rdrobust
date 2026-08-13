## Frozen numerical baselines.
##
## These pin the exact output of the current implementation so that any change
## in the numbers is deliberate rather than accidental. A failure here is not
## automatically a bug: it means the estimator moved, and the move must be
## explained (and the baseline updated in the same commit as the fix).
##
## The RDsenate baselines use the shipped dataset and therefore carry no
## dependence on the random number generator or the R version.

test_that("RDsenate fixed-bandwidth baseline is stable", {
  data(rdrobust_RDsenate, envir = environment())
  y <- rdrobust_RDsenate$vote
  x <- rdrobust_RDsenate$margin

  est <- rdrobust(y = y, x = x, h = 15)

  expect_equal(as.numeric(est$coef),
               c(7.48728585809499, 9.08562818492005, 9.08562818492005),
               tolerance = 1e-10)
  expect_equal(as.numeric(est$se),
               c(1.55973178848191, 1.55973178848191, 2.24067206402775),
               tolerance = 1e-10)
  expect_equal(as.numeric(est$N_h), c(319, 288))
  expect_equal(as.numeric(est$N), c(595, 702))
})

test_that("RDsenate data-driven baseline is stable", {
  data(rdrobust_RDsenate, envir = environment())
  y <- rdrobust_RDsenate$vote
  x <- rdrobust_RDsenate$margin

  est <- rdrobust(y = y, x = x)

  expect_equal(as.numeric(est$bws),
               c(17.7543972211599, 28.0280871511408,
                 17.7543972211599, 28.0280871511408),
               tolerance = 1e-9)
  expect_equal(as.numeric(est$coef),
               c(7.41413080197352, 7.50650247058758, 7.50650247058758),
               tolerance = 1e-9)
  expect_equal(as.numeric(est$se),
               c(1.45871602382597, 1.45871602382597, 1.74125841460802),
               tolerance = 1e-9)
  expect_equal(as.numeric(est$N_h), c(360, 323))
})

test_that("RDsenate bandwidth-selector baselines are stable", {
  data(rdrobust_RDsenate, envir = environment())
  y <- rdrobust_RDsenate$vote
  x <- rdrobust_RDsenate$margin

  bw <- rdbwselect(y = y, x = x, all = TRUE)

  expect_equal(as.numeric(bw$bws["mserd", ]),
               c(17.7543972211599, 17.7543972211599,
                 28.0280871511408, 28.0280871511408),
               tolerance = 1e-9)
  expect_equal(as.numeric(bw$bws["msetwo", ]),
               c(16.1698202033388, 18.1264623349676,
                 27.1038901767765, 29.3435512423602),
               tolerance = 1e-9)
  expect_equal(as.numeric(bw$bws["cerrd", ]),
               c(12.4067757744621, 12.4067757744621,
                 28.0280871511408, 28.0280871511408),
               tolerance = 1e-9)
})

test_that("RDsenate kink baseline with a non-default kernel and vce is stable", {
  data(rdrobust_RDsenate, envir = environment())
  y <- rdrobust_RDsenate$vote
  x <- rdrobust_RDsenate$margin

  est <- rdrobust(y = y, x = x, h = 15, p = 2, deriv = 1,
                  kernel = "epa", vce = "hc3")

  expect_equal(as.numeric(est$coef),
               c(1.91767244605116, 3.31534165956668, 3.31534165956668),
               tolerance = 1e-9)
  expect_equal(as.numeric(est$se),
               c(0.779215462267211, 0.779215462267211, 1.86452219485337),
               tolerance = 1e-9)
})

test_that("covariate-adjusted clustered baseline is stable", {
  fx <- make_rd_fixture()
  est <- rdrobust(y = fx$y, x = fx$x, covs = fx$Z, cluster = fx$cluster,
                  vce = "cr1", h = 0.3)

  expect_equal(as.numeric(est$coef),
               c(0.588602499100754, 0.563768246801708, 0.563768246801708),
               tolerance = 1e-9)
  expect_equal(as.numeric(est$se),
               c(0.222366879403164, 0.222366879403164, 0.296682179015828),
               tolerance = 1e-9)
})

test_that("fuzzy baseline is stable", {
  fz <- make_fuzzy_fixture()
  est <- rdrobust(y = fz$y, x = fz$x, fuzzy = fz$d, h = 0.3)

  expect_equal(as.numeric(est$coef),
               c(1.15411192106143, 1.30906992078803, 1.30906992078803),
               tolerance = 1e-9)
  expect_equal(as.numeric(est$se),
               c(0.254478140749401, 0.254478140749401, 0.353555726921529),
               tolerance = 1e-9)
})

test_that("RDsenate rdplot bin counts are stable", {
  data(rdrobust_RDsenate, envir = environment())
  pl <- rdplot(y = rdrobust_RDsenate$vote, x = rdrobust_RDsenate$margin,
               hide = TRUE)
  expect_equal(as.numeric(pl$J), c(15, 35))
})
