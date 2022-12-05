\name{rdrobust-package}
\alias{rdrobust-package}
\docType{package}
\title{Robust Data-Driven Statistical Inference in RD Designs}
\description{
   Regression-discontinuity (RD) designs are quasi-experimental research designs popular in social, behavioral and natural sciences. The RD design is usually employed to study the (local) causal effect of a treatment, intervention or policy. This package provides tools for data-driven graphical and analytical statistical inference in RD designs: \code{\link{rdrobust}} to construct local-polynomial point estimators and robust confidence intervals for average treatment effects at the cutoff in Sharp, Fuzzy and Kink RD settings, \code{\link{rdbwselect}} to perform bandwidth selection for the different procedures implemented, and \code{\link{rdplot}} to conduct exploratory data analysis (RD plots).
}
\details{
\tabular{ll}{
Package: \tab rdrobust\cr
Type: \tab Package\cr
Version: \tab 2.1.1\cr 
Date: \tab 2022-12-04\cr
License: \tab GPL-2\cr
}
Function for statistical inference: \code{\link{rdrobust}}\cr
Function for bandwidths selection: \code{\link{rdbwselect}}\cr
Function for exploratory data analysis (RD plots): \code{\link{rdplot}}
}
\author{
Sebastian Calonico, Columbia University, New York, NY.  \email{sebastian.calonico@columbia.edu}.

Matias D. Cattaneo, Princeton University, Princeton, NJ.  \email{cattaneo@princeton.edu}.

Max H. Farrell, University of Chicago, Chicago, IL. \email{max.farrell@chicagobooth.edu}.

Rocio Titiunik, Princeton University, Princeton, NJ.  \email{titiunik@princeton.edu}.
}
