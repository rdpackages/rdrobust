********************************************************************************
* RDROBUST STATA PACKAGE -- rdrobustplot
* Diagnostic plot for a previous rdrobust result. Mirrors the plot.rdrobust()
* S3 method in the R package and plot_rdrobust() in the Python package.
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell, Rocio Titiunik
********************************************************************************
*! version 11.1.0 22may2026

capture program drop rdrobustplot
program define rdrobustplot, rclass
	version 16.0
	* ST-9: scale/xlabel/ylabel/xtitle/ytitle were accepted and then never
	* referenced -- silently dead. They are ordinary twoway options, so they
	* are now folded into the graph_options() handed to rdplot (see below).
	* col_dots()/col_lines() are NOT here any more: they are per-plot marker
	* and line colours, and rdplot builds its own twoway call with no hook for
	* them, so they could never have worked. Removed rather than left as
	* silent no-ops -- Stata now reports "option col_dots() not allowed".
	syntax [, nbins(string) binselect(string) NOCI scale(string) ///
		title(string asis) xlabel(string) ylabel(string) ///
		xtitle(string asis) ytitle(string asis) ///
		shade * ]

	* -----------------------------------------------------------------------
	* Validate: require a previous rdrobust call
	* -----------------------------------------------------------------------
	if ("`e(cmd)'" != "rdrobust") {
		di as error "{err}{cmd:rdrobustplot} must follow a {cmd:rdrobust} call"
		exit 301
	}

	* -----------------------------------------------------------------------
	* Pull relevant quantities from e(). Note: `local x = e(...)` with the `=`
	* sign evaluates missing-ereturn-locals to ".", so for strings that might
	* be unset (covs, kernel) we use macro expansion instead.
	* -----------------------------------------------------------------------
	local y       "`e(depvar)'"
	local x       "`e(runningvar)'"
	local covs    "`e(covs)'"
	local c       = e(c)
	local p       = e(p)
	local h_l     = e(h_l)
	local h_r     = e(h_r)
	local tau     = e(tau_cl)
	local se_rb   = e(se_tau_rb)
	local ci_l    = e(ci_l_rb)
	local ci_r    = e(ci_r_rb)
	local level   = e(level)
	local kernel  = lower(substr("`e(kernel)'", 1, 3))   /* "Triangular" -> "tri" */

	* Guard against missing e() values that would otherwise propagate NaN
	* through the subtitle's p-value computation.
	if (missing(`tau') | missing(`se_rb') | `se_rb' <= 0) {
		di as error "{err}{cmd:rdrobustplot} cannot build subtitle: e(tau_cl) or e(se_tau_rb) is missing or non-positive"
		exit 198
	}

	* -----------------------------------------------------------------------
	* Defaults
	* -----------------------------------------------------------------------
	if ("`nbins'"   == "") local nbins     = "20 20"
	if ("`binselect'"=="") local binselect = "esmv"
	if ("`scale'"   == "") local scale     = ""

	* Build effect annotation (shown as subtitle above the plot)
	local pv = 2*normal(-abs(`tau'/`se_rb'))
	local stars = ""
	if (`pv' < 0.10) local stars = "*"
	if (`pv' < 0.05) local stars = "**"
	if (`pv' < 0.01) local stars = "***"

	local ci_lo_s = strofreal(`ci_l', "%7.3f")
	local ci_hi_s = strofreal(`ci_r', "%7.3f")
	local tau_s   = strofreal(`tau',  "%7.3f")
	local lev_s   = strofreal(`level',"%2.0f")

	local subtitle = "RD = `tau_s'`stars'  (`lev_s'% RBC CI: [`ci_lo_s', `ci_hi_s'])"

	if (`"`title'"' == "") local title `"RD plot: `y' vs `x'"'

	* -----------------------------------------------------------------------
	* Delegate to rdplot using the SAME bandwidth the rdrobust call used
	* -----------------------------------------------------------------------
	local ci_flag = cond("`noci'"!="", "", "ci(`level')")
	local shade_flag = cond("`shade'"!="", "shade", "")
	local covs_flag  = cond("`covs'"!="",  "covs(`covs')", "")

	* ST-9: build the twoway option list so the previously-dead cosmetic
	* options actually reach the graph.
	local gopts `"title(`"`title'"') subtitle(`"`subtitle'"')"'
	if (`"`xtitle'"' != "") local gopts `"`gopts' xtitle(`"`xtitle'"')"'
	if (`"`ytitle'"' != "") local gopts `"`gopts' ytitle(`"`ytitle'"')"'
	if ("`xlabel'"  != "")  local gopts `"`gopts' xlabel(`xlabel')"'
	if ("`ylabel'"  != "")  local gopts `"`gopts' ylabel(`ylabel')"'
	if ("`scale'"   != "")  local gopts `"`gopts' scale(`scale')"'

	* ST-9: rdplot is eclass, so it CLEARS e() -- after one rdrobustplot the
	* e(cmd)=="rdrobust" guard at the top of this program failed and a SECOND
	* call died rc=301. Park the rdrobust estimates across the delegation and
	* put them back, so rdrobustplot is repeatable and leaves e() as it found
	* it.
	tempname _rdrp_est
	_estimates hold `_rdrp_est', nullok
	capture noisily rdplot `y' `x' , c(`c') h(`h_l' `h_r') p(`p') ///
		nbins(`nbins') binselect(`binselect') kernel(`kernel') ///
		`covs_flag' `ci_flag' `shade_flag' ///
		graph_options(`gopts') ///
		`options'
	local _rdp_rc = _rc
	_estimates unhold `_rdrp_est'
	if (`_rdp_rc') exit `_rdp_rc'

	* -----------------------------------------------------------------------
	* Return the annotation bits so downstream scripts can reuse
	* -----------------------------------------------------------------------
	return local subtitle = `"`subtitle'"'
	return scalar tau    = `tau'
	return scalar se_rb  = `se_rb'
	return scalar ci_l   = `ci_l'
	return scalar ci_r   = `ci_r'
	return scalar pvalue = `pv'
end
