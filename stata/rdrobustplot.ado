********************************************************************************
* RDROBUST STATA PACKAGE -- rdrobustplot
* Diagnostic plot for a previous rdrobust result. Mirrors the plot.rdrobust()
* S3 method in the R package and plot_rdrobust() in the Python package.
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell, Rocio Titiunik
********************************************************************************
*! version 11.1.1 01oct2026

capture program drop rdrobustplot
program define rdrobustplot, rclass
	version 16.0
	* ST-9: scale/xlabel/ylabel/xtitle/ytitle were accepted and then never
	* referenced -- silently dead. They are ordinary twoway options, so they
	* are now folded into the graph_options() handed to rdplot (see below).
	* col_dots()/col_lines() are per-plot marker and line colours, and rdplot
	* builds its own twoway call with no hook for them, so they have never had
	* an effect. They are still accepted, so existing scripts keep running,
	* with a note saying they are ignored.
	* Two passthrough routes, kept separate on purpose:
	*  - rdplot ANALYSIS options (masspoints() etc.) are declared by name
	*    below and forwarded to rdplot as such;
	*  - the trailing `*' collects everything else (legend(), xline(), ...),
	*    which is treated as twoway options and folded into graph_options().
	* The wildcard used to be forwarded to rdplot as top-level options,
	* which rdplot silently swallowed before ST-10 removed its own wildcard
	* and rejects with rc=198 after -- so there was never a working route
	* for twoway options until now. graph_options() is declared explicitly
	* so a caller using rdplot's spelling merges instead of nesting.
	syntax [, nbins(string) binselect(string) NOCI scale(string) ///
		title(string asis) xlabel(string) ylabel(string) ///
		xtitle(string asis) ytitle(string asis) ///
		masspoints(string) covs_drop(string) covs_eval(string) ///
		support(string) genvars nochecks PRECision(string) ///
		graph_options(string asis) col_dots(string) col_lines(string) shade * ]

	if (`"`col_dots'`col_lines'"' != "") di as text "Note: col_dots() and col_lines() have no effect and are ignored."

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
	* NOCI is declared in capitals (no abbreviation), so its macro is `noci'.
	* The lowercase nochecks below is a negated option, stored in `checks'.
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
	* Route explicit graph_options() and any wildcard-collected twoway
	* options into the same slot. rdplot appends graph_options() AFTER its
	* own twoway settings, so later-wins lets e.g. legend(off) override the
	* default legend. A typo'd rdplot option lands here too and errors from
	* -twoway- ("option ... not allowed") rather than from rdplot.
	if (`"`graph_options'"' != "") local gopts `"`gopts' `graph_options'"'
	if (`"`options'"'       != "") local gopts `"`gopts' `options'"'

	* ST-9: rdplot is eclass, so it CLEARS e() -- after one rdrobustplot the
	* e(cmd)=="rdrobust" guard at the top of this program failed and a SECOND
	* call died rc=301. Park the rdrobust estimates across the delegation and
	* put them back, so rdrobustplot is repeatable and leaves e() as it found
	* it.
	tempname _rdrp_est
	_estimates hold `_rdrp_est', nullok
	* Analysis options forwarded to rdplot by name (empty ones are inert).
	local rdp_pass ""
	if ("`masspoints'" != "") local rdp_pass "`rdp_pass' masspoints(`masspoints')"
	if ("`covs_drop'"  != "") local rdp_pass "`rdp_pass' covs_drop(`covs_drop')"
	if ("`covs_eval'"  != "") local rdp_pass "`rdp_pass' covs_eval(`covs_eval')"
	if ("`support'"    != "") local rdp_pass "`rdp_pass' support(`support')"
	if ("`precision'"  != "") local rdp_pass "`rdp_pass' precision(`precision')"
	capture noisily rdplot `y' `x' , c(`c') h(`h_l' `h_r') p(`p') ///
		nbins(`nbins') binselect(`binselect') kernel(`kernel') ///
		`covs_flag' `ci_flag' `shade_flag' `genvars' `checks' `rdp_pass' ///
		graph_options(`gopts')
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
