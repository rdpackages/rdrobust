# Changelog

## October 1, 2026 update: R 4.1.1, Python 2.1.1, Stata 11.1.1

Fixes from a three-language audit. R `4.1.1` replaces the unreleased GitHub
`4.1.0` as the next CRAN version; Python `2.1.1` and Stata `11.1.1` carry the
same fixes.

Changes in results:

- Sharp RD with covariates and `deriv >= 2`: apply the missing
  `factorial(deriv)` to the point estimate (R, Python, Stata). With
  `deriv = 2` the estimate was half its correct value.
- Python: `deriv >= 2` no longer crashes.
- `rdbwselect()` with `stdvars = TRUE` (the default) reports the effective
  sample size and the cutoff on the original scale (R, Python).
- R: fuzzy designs with `sharpbw = TRUE` and a nonzero cutoff no longer crash.
- `rdplot()`: bin edges stay with their own bins when a side has empty bins
  (R, Python); Stata `genvars` describe the estimation sample, and Stata's
  per-bin confidence intervals use `N - 1` degrees of freedom, as R and Python.
- Bandwidth selection: the `bwcheck` floor and the first pilot bandwidth are
  padded by a factor `1 + sqrt(eps)`, so the observation that sets them keeps a
  positive kernel weight (R, Python, Stata). Selected bandwidths move by less
  than `1e-6` in relative terms.
- R: the default `ginv.tol` is `1e-15` (was `1e-20`), as in Python and Stata.
  With a covariate that is nearly collinear on one side of the cutoff, the
  selected bandwidth no longer depends on row order or numerical noise.

Errors and warnings:

- When the bandwidth-selection pilots are not identified (fewer distinct values
  of `x` in a pilot window than polynomial coefficients, as with a running
  variable that takes a few values), all three languages stop with "Not enough
  variability". Before, R returned a bandwidth from a generalized inverse, NaN
  or an unrelated error; Python raised a linear-algebra error; and Stata printed
  a message and continued, usually with the full range as bandwidth.
- A side of the cutoff with fewer than `q + 1` distinct values of `x`
  (`p + 1` in `rdbwselect()`) stops with an informative error.
- With clustered standard errors, a warning when either side has fewer than 10
  clusters within the bandwidth, and when a standard error is exactly 0.
- Input validation: `vce`, `kernel` and `bwselect` are case-insensitive in
  `rdbwselect()` as in `rdrobust()`; `nnmatch`, negative weights,
  `binselect`, collinear covariates with `covs_drop = FALSE`, and a cutoff
  outside the range of `x` give informative errors.
- Stata: `e(sample)` is set, user matrices and tempvars are no longer
  overwritten, and `rdrobustplot` accepts twoway options.

Tests: R has a `testthat` suite (116 tests); Python has 37 tests.

## September 27, 2026 update

Python `2.1.0` is published on PyPI. The corresponding R fixes are available
in GitHub version `4.1.0`; CRAN submission is pending. Stata fixes are available
from GitHub.

- Python `2.1.0` also includes the earlier changes enabling standardization by
  default during automatic bandwidth selection and improving scale stability
  for nearest-neighbor ties and mass-point bandwidth floors.
- Respect supplied estimation and bias bandwidths for samples with fewer than
  20 observations, and honor `rho` consistently in R, Python, and Stata.
- Align automatic small-sample bandwidths across implementations using the
  maximum distance from the cutoff. Fix Python's scalar-maximum error and
  preserve original running-variable units in R when standardization is enabled.
- Reject local polynomial fits without enough distinct, positively weighted
  running-variable values on each side. R now stops when polynomial inversion
  fails, preventing misleading inference from a generalized inverse.
- Accept Python list and tuple bandwidths with `rho`, validate bandwidth inputs,
  and correctly apply integer-valued floating-point subset indices.

## Modernization Summary: May 11-17, 2026

- Prepared the `rdrobust` 4.0 release materials across the repository, including refreshed package metadata, public README content, examples, and bundled replication data.
- Modernized the Python package by cleaning generated build and editor artifacts, refreshing package metadata, bundling the Senate dataset with a loader, adding `plot_rdrobust`, and adding focused numerical contract tests.
- Modernized the R package metadata, namespace, documentation, bundled data, and package build configuration, including support files for the 4.0 release.
- Refreshed Stata release artifacts, including ado files, help files, PDFs, package metadata, compiled Mata modules, the Senate dataset, updated illustrations, and the `rdrobustplot` postestimation diagnostic plot.
- Added GitHub-facing project infrastructure: CI for repository layout, Python package checks, and R package checks; Python publishing workflow; Dependabot configuration; issue and pull request templates; and a security reporting policy.
