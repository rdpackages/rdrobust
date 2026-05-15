\name{plot.rdrobust}
\alias{plot.rdrobust}
\title{Plot Method for rdrobust Objects}
\description{
Produces a visual summary of an \code{\link{rdrobust}} result. The main panel
shows a binned scatter plot of the outcome against the running variable within
the estimated bandwidth, overlaid with the local-polynomial fit curves and
(optionally) pointwise confidence bands. Bin sizes are scaled to the average
kernel weight of the observations they contain, so bins closest to the cutoff
appear largest. An optional effect panel below the main plot displays the
conventional point estimate with its robust bias-corrected confidence interval
and significance stars.
}
\usage{
\method{plot}{rdrobust}(x, y, x_run, nbins = 20,
     ci = TRUE, show_effect = FALSE,
     title = NULL,
     x.label = "Running Variable",
     y.label = "Outcome",
     x.lim = NULL, y.lim = NULL,
     col.l = "#3B7DD8", col.r = "#D95F3B",
     base_size = 14, \dots)
}
\arguments{
  \item{x}{an object of class \code{"rdrobust"}, as returned by \code{\link{rdrobust}}.}
  \item{y}{numeric vector of outcome values (same length and order as used in the original \code{rdrobust} call, before any \code{subset} was applied).}
  \item{x_run}{numeric vector of running variable values (same length and order as \code{y}).}
  \item{nbins}{number of bins per side used for the binned scatter plot. Default is \code{20}.}
  \item{ci}{logical. If \code{TRUE} (default), shaded pointwise confidence bands are drawn around each polynomial fit curve. The bands use the conventional variance-covariance matrices \code{V_cl_l} and \code{V_cl_r} stored in the \code{rdrobust} object, evaluated at the confidence level stored in the \code{rdrobust} object.}
  \item{show_effect}{logical. If \code{TRUE}, a second panel is printed below the main plot showing the conventional point estimate with its robust bias-corrected (RBC) confidence interval as a horizontal error bar, together with significance stars based on the robust p-value. Combining both panels requires the \pkg{gridExtra} package; if unavailable the two plots are printed sequentially. Default is \code{FALSE}.}
  \item{title}{character string for the plot title. Default is \code{NULL} (no title).}
  \item{x.label}{label for the horizontal axis. Default is \code{"Running Variable"}.}
  \item{y.label}{label for the vertical axis. Default is \code{"Outcome"}.}
  \item{x.lim}{numeric vector of length 2 specifying the horizontal axis limits. Default is \code{NULL} (automatic).}
  \item{y.lim}{numeric vector of length 2 specifying the vertical axis limits. Default is \code{NULL} (automatic).}
  \item{col.l}{color for the left-side (control) elements. Default is \code{"#3B7DD8"}.}
  \item{col.r}{color for the right-side (treated) elements. Default is \code{"#D95F3B"}.}
  \item{base_size}{base font size in points passed to the underlying \pkg{ggplot2} theme. All text elements scale relative to this value. Default is \code{14}.}
  \item{\dots}{additional arguments (currently unused).}
}
\details{
Only observations within the main estimation bandwidth \code{[c - h_l, c + h_r]}
are shown. Each side is divided into \code{nbins} equal-width bins; the plotted
point for each bin is the bin mean of both \code{x_run} and \code{y}. The size
of each bin's dot is proportional to the average kernel weight of the
observations it contains, computed from the kernel stored in
\code{x$kernel} (triangular, Epanechnikov, or uniform) evaluated at the
normalized distance \code{(x_i - c) / h}.

The polynomial fit curve on each side is the fitted conditional mean implied by
the local-polynomial coefficients \code{beta_Y_p_l} and \code{beta_Y_p_r}. The
pointwise confidence bands (\code{ci = TRUE}) are computed as
\deqn{\hat{\mu}(t) \pm z_{\alpha/2} \sqrt{\mathbf{x}(t)^\top V \mathbf{x}(t)}}
where \eqn{\mathbf{x}(t) = (1, t-c, \ldots, (t-c)^p)^\top} and \eqn{V} is
the conventional variance-covariance matrix of the polynomial coefficients.

When \code{show_effect = TRUE} the effect panel reports:
\itemize{
  \item Point estimate: conventional (\code{coef["Conventional"]}).
  \item Interval: robust bias-corrected (\code{ci["Robust",]}).
  \item Stars: based on \code{pv["Robust",]}: \code{*} p<0.10, \code{**} p<0.05, \code{***} p<0.01.
}
}
\value{
Invisibly returns the \pkg{ggplot2} plot object when \code{show_effect = FALSE},
or a named list \code{list(rd_plot, effect_plot)} of two \pkg{ggplot2} objects
when \code{show_effect = TRUE}. In both cases the plot(s) are printed as a
side-effect.
}
\references{
Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014.
\href{https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf}{Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs.}
\emph{Econometrica} 82(6): 2295-2326.

Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2019.
\href{https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf}{Regression Discontinuity Designs using Covariates.}
\emph{Review of Economics and Statistics}, 101(3): 442-451.
}
\author{
Sebastian Calonico, University of California, Davis, CA.  \email{scalonico@ucdavis.edu}.

Matias D. Cattaneo, Princeton University, Princeton, NJ.  \email{matias.d.cattaneo@gmail.com}.

Max H. Farrell, University of California, Santa Barbara, CA. \email{mhfarrell@gmail.com}.

Rocio Titiunik, Princeton University, Princeton, NJ.  \email{rocio.titiunik@gmail.com}.
}
\examples{
x <- runif(500, -1, 1)
y <- 5 + 3*x + 2*(x >= 0) + rnorm(500)
r <- rdrobust(y, x)

## Main RD plot with CI bands
plot(r, y, x)

## Add effect panel
plot(r, y, x, show_effect = TRUE)
}
\seealso{
  \code{\link{rdrobust}}, \code{\link{rdbwselect}}, \code{\link{rdplot}}
}
\keyword{RDD}
\keyword{Plots}
