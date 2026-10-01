# Changelog

## October 1, 2026 update

Fixes from a three-language audit, added to R `4.1.0` before its CRAN
submission. Python `2.1.0` on PyPI does not include them yet.

- Sharp RD with covariates and `deriv >= 2`: apply the missing
  `factorial(deriv)` to the point estimate (R, Python, Stata). With
  `deriv = 2` the estimate was half its correct value.
- Python: `deriv >= 2` no longer crashes.
- `rdbwselect()` with `stdvars = TRUE` (the default) reports the effective
  sample size and the cutoff on the original scale (R, Python).
- R: fuzzy designs with `sharpbw = TRUE` and a nonzero cutoff no longer crash.
- `rdplot()`: bin edges stay with their own bins when a side has empty bins
  (R, Python); Stata `genvars` describe the estimation sample.
- Input validation: `vce`, `kernel` and `bwselect` are case-insensitive in
  `rdbwselect()` as in `rdrobust()`; `nnmatch`, negative weights,
  `binselect`, collinear covariates with `covs_drop = FALSE`, and a cutoff
  outside the range of `x` give informative errors.
- Stata: `e(sample)` is set, user matrices and tempvars are no longer
  overwritten, and `rdrobustplot` accepts twoway options.
- R: added a `testthat` suite (111 tests).

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
