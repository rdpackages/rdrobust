{smcl}
{* *!version 11.1.0  2026-05-22}{...}
{viewerjumpto "Syntax" "rdrobustplot##syntax"}{...}
{viewerjumpto "Description" "rdrobustplot##description"}{...}
{viewerjumpto "Options" "rdrobustplot##options"}{...}
{viewerjumpto "Stored results" "rdrobustplot##results"}{...}
{viewerjumpto "Examples" "rdrobustplot##examples"}{...}

{title:Title}

{p 4 8}{cmd:rdrobustplot} {hline 2} Diagnostic plot for a previous {cmd:rdrobust} result.{p_end}

{marker syntax}{...}
{title:Syntax}

{p 4 8}{cmd:rdrobustplot}
[{cmd:,}
{cmd:nbins(}{it:# #}{cmd:)}
{cmd:binselect(}{it:binmethod}{cmd:)}
{cmd:noci}
{cmd:shade}
{cmd:scale(}{it:# #}{cmd:)}
{cmd:title(}{it:string}{cmd:)}
{cmd:xtitle(}{it:string}{cmd:)}
{cmd:ytitle(}{it:string}{cmd:)}
{cmd:xlabel(}{it:rule}{cmd:)}
{cmd:ylabel(}{it:rule}{cmd:)}
{it:rdplot_options}
]{p_end}

{marker description}{...}
{title:Description}

{p 4 8}{cmd:rdrobustplot} produces a diagnostic RD plot for the results of the
most recent {help rdrobust:rdrobust} call. It wraps {help rdplot:rdplot} with the
main bandwidth and polynomial order that {cmd:rdrobust} used, and decorates the
subtitle with the estimated RD coefficient and its robust
bias-corrected confidence interval.{p_end}

{p 4 8}This mirrors the {cmd:plot.rdrobust()} S3 method in the R package and the
{cmd:plot_rdrobust()} function in the Python package.{p_end}

{p 4 8}{it:Requires Stata 16 or later.}{p_end}

{marker options}{...}
{title:Options}

{p 4 8}{cmd:nbins(}{it:# #}{cmd:)} number of bins per side; default is {cmd:20 20}.{p_end}

{p 4 8}{cmd:binselect(}{it:binmethod}{cmd:)} bin selection rule; default is
{cmd:esmv} (mimicking-variance evenly-spaced bins). See {help rdplot}.{p_end}

{p 4 8}{cmd:noci} suppresses the per-bin confidence intervals.{p_end}

{p 4 8}{cmd:shade} draws the pointwise confidence bands as a shaded ribbon
instead of error bars.{p_end}

{p 4 8}{cmd:title()}, {cmd:xtitle()}, {cmd:ytitle()}, {cmd:xlabel()},
{cmd:ylabel()}, {cmd:scale()} are forwarded to the underlying graph. The
default title follows {cmd:rdrobust}'s outcome and running variable.{p_end}

{p 4 8}{cmd:masspoints()}, {cmd:covs_drop()}, {cmd:covs_eval()},
{cmd:support()}, {cmd:genvars}, {cmd:nochecks} and {cmd:precision()} are
forwarded to {help rdplot:rdplot} as analysis options.{p_end}

{p 4 8}{cmd:graph_options(}{it:twoway options}{cmd:)}, and any option not
listed above, are appended to the underlying {help twoway} call after the
defaults, so e.g. {cmd:legend(off)} or {cmd:xline(0.5)} take effect
(later-wins). A misspelled option therefore errors from the graph command
("option ... not allowed") rather than silently disappearing.{p_end}

{p 4 8}{it:Note:} {cmd:col_dots()} and {cmd:col_lines()} were accepted by
earlier versions but never had any effect -- the binned means and the fit line
are drawn inside {cmd:rdplot}, which exposes no hook for per-plot colours.
They are no longer accepted, so the mistake is reported instead of being
silently ignored.{p_end}

{marker results}{...}
{title:Stored results}

{p 4 8}{cmd:rdrobustplot} is {cmd:r}-class and stores the annotation shown in
the subtitle, so a script can reuse the same numbers:{p_end}

{synoptset 20 tabbed}{...}
{p2col 5 20 24 2: Scalars}{p_end}
{synopt:{cmd:r(tau)}}RD point estimate, {cmd:e(tau_cl)} of the preceding {cmd:rdrobust}{p_end}
{synopt:{cmd:r(se_rb)}}robust bias-corrected standard error{p_end}
{synopt:{cmd:r(ci_l)}}lower robust bias-corrected confidence limit{p_end}
{synopt:{cmd:r(ci_r)}}upper robust bias-corrected confidence limit{p_end}
{synopt:{cmd:r(pvalue)}}two-sided p-value implied by {cmd:r(tau)} and {cmd:r(se_rb)}{p_end}

{p2col 5 20 24 2: Macros}{p_end}
{synopt:{cmd:r(subtitle)}}the annotation string drawn above the plot{p_end}

{p 4 8}The {cmd:e()} results of the preceding {cmd:rdrobust} call are left
intact, so {cmd:rdrobustplot} may be called repeatedly (and followed by other
post-estimation commands). Earlier versions delegated to the {cmd:e}-class
{cmd:rdplot}, which cleared {cmd:e()} and made a second {cmd:rdrobustplot}
call fail with r(301).{p_end}

{marker examples}{...}
{title:Examples}

{phang2}{stata sysuse rdrobust_senate, clear}{p_end}
{phang2}{stata rdrobust vote margin}{p_end}
{phang2}{stata rdrobustplot}{p_end}

{pstd}With cluster-robust variance and fuzzy RD:{p_end}
{phang2}{stata rdrobust vote margin, vce(cr3 state)}{p_end}
{phang2}{stata rdrobustplot, shade}{p_end}

{title:Authors}

{p 4 8}Sebastian Calonico, University of California, Davis. {browse "mailto:scalonico@ucdavis.edu":scalonico@ucdavis.edu}.{p_end}
{p 4 8}Matias D. Cattaneo, Princeton University. {browse "mailto:matias.d.cattaneo@gmail.com":matias.d.cattaneo@gmail.com}.{p_end}
{p 4 8}Max H. Farrell, University of California, Santa Barbara. {browse "mailto:mhfarrell@gmail.com":mhfarrell@gmail.com}.{p_end}
{p 4 8}Rocio Titiunik, Princeton University. {browse "mailto:rocio.titiunik@gmail.com":rocio.titiunik@gmail.com}.{p_end}
