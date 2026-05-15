{smcl}
{* *!version 11.0.0  2026-05-13}{...}
{viewerjumpto "Syntax" "rdrobustplot##syntax"}{...}
{viewerjumpto "Description" "rdrobustplot##description"}{...}
{viewerjumpto "Options" "rdrobustplot##options"}{...}
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
{cmd:col_dots(}{it:color}{cmd:)}
{cmd:col_lines(}{it:color}{cmd:)}
{it:graph_options}
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

{p 4 8}{cmd:title()}, {cmd:xtitle()}, {cmd:ytitle()} custom text for the plot
elements (defaults follow {cmd:rdrobust}'s outcome / running variable).{p_end}

{p 4 8}{cmd:col_dots()}, {cmd:col_lines()} colors for the binned means and the
polynomial fit line.{p_end}

{p 4 8}Any other option is passed through to {cmd:rdplot}'s
{cmd:graph_options()}.{p_end}

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
