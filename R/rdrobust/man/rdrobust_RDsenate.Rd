\name{rdrobust_RDsenate}
\alias{rdrobust_RDsenate}
\docType{data}
\title{RD Senate Data}
\description{
Extract of the dataset constructed by Cattaneo, Frandsen, and Titiunik (2015), which includes measures of incumbency advantage in the U.S. Senate for the period 1914-2010. The dataset contains the running variable (Democratic vote margin at election \eqn{t}), the main outcome (Democratic vote share at election \eqn{t+2}), additional outcomes and lagged vote shares useful as covariates, and a state identifier suitable for cluster-robust inference.
}
\usage{data(rdrobust_RDsenate)}
\format{
  A data frame with 1390 observations on the following 17 variables.
  \describe{
    \item{\code{state}}{character. U.S. state name; can be used as a cluster variable for cluster-robust variance estimation.}
    \item{\code{year}}{numeric. Election year.}
    \item{\code{margin}}{numeric. Democratic vote margin at election \eqn{t} (running variable); equals \code{demmv} in the original dataset.}
    \item{\code{vote}}{numeric. Democratic vote share at election \eqn{t+2} (main outcome); equals \code{demvoteshfor2} in the original dataset.}
    \item{\code{class}}{numeric. U.S. Senate class (1, 2, or 3).}
    \item{\code{termshouse}}{numeric. Number of terms served in the U.S. House.}
    \item{\code{termssenate}}{numeric. Number of terms served in the U.S. Senate.}
    \item{\code{dopen}}{numeric. Indicator for open-seat election.}
    \item{\code{population}}{numeric. State population.}
    \item{\code{presdemvoteshlag1}}{numeric. Lagged Democratic presidential vote share.}
    \item{\code{demvoteshlag1}}{numeric. Democratic vote share at election \eqn{t} (same cycle as running variable); useful as a covariate.}
    \item{\code{demvoteshlag2}}{numeric. Democratic vote share at election \eqn{t-2}; useful as a covariate.}
    \item{\code{demvoteshfor1}}{numeric. Democratic vote share at election \eqn{t+1}.}
    \item{\code{demwinprv1}}{numeric. Indicator for Democratic win at election \eqn{t-1}.}
    \item{\code{demwinprv2}}{numeric. Indicator for Democratic win at election \eqn{t-2}.}
    \item{\code{dmidterm}}{numeric. Indicator for midterm election year.}
    \item{\code{dpresdem}}{numeric. Indicator for Democratic president in office.}
  }
}
\source{
    Cattaneo, M. D., Frandsen, B., and R. Titiunik. 2015. \href{https://rdpackages.github.io/references/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate.} \emph{Journal of Causal Inference} 3(1): 1-24.
}
\references{
    Cattaneo, M. D., Frandsen, B., and R. Titiunik. 2015. \href{https://rdpackages.github.io/references/Cattaneo-Frandsen-Titiunik_2015_JCI.pdf}{Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate.} \emph{Journal of Causal Inference} 3(1): 1-24.
}
\keyword{datasets}
