********************************************************************************
** RDROBUST Stata Package
** Illustration of new features (v11.0.0, 2026-05-13)
********************************************************************************
** New in this release:
**   (1) Full HC family of heteroskedastic variance estimators:
**         vce(hc0) vce(hc1) vce(hc2) vce(hc3)
**   (2) Cluster-robust variance estimators CR1 / CR2 / CR3:
**         vce(cr1 clustervar)  -- df correction (equivalent to vce(cluster ...))
**         vce(cr2 clustervar)  -- Bell-McCaffrey (2002)
**         vce(cr3 clustervar)  -- Pustejovsky-Tipton (2018), approximately unbiased
**   (3) rdrobustplot -- diagnostic plot for a previous rdrobust result
********************************************************************************
clear all
set more off
set linesize 100

********************************************************************************
** Setup: load the unified senate dataset (bundled with the rdrobust SSC package)
********************************************************************************
findfile rdrobust_senate.dta
use "`r(fn)'", clear
describe, short

********************************************************************************
** (1) HC family of heteroskedastic variance estimators
********************************************************************************
di _n "======================================================================"
di       "(1) Heteroskedastic-robust variance: NN vs HC0-HC3"
di       "======================================================================"

foreach v in nn hc0 hc1 hc2 hc3 {
    qui rdrobust vote margin, vce(`v')
    display as text "vce=" %-5s "`v'" "   tau_cl = " %7.4f e(tau_cl) ///
                    "   se_rb = " %7.4f e(se_tau_rb) ///
                    "   pv_rb = " %6.4f e(pv_rb)
}

********************************************************************************
** (2) Cluster-robust variance: CR1 / CR2 / CR3
********************************************************************************
di _n "======================================================================"
di       "(2) Cluster-robust variance: CR family (clustered on state)"
di       "======================================================================"
di _n "vce(cluster state) is an alias for vce(cr1 state). CR2 is Bell-McCaffrey;"
di       "CR3 is Pustejovsky-Tipton (approximately unbiased with few clusters)."

* Note: `state' is a string variable; vce(cluster ...) handles it directly.
foreach v in cr1 cr2 cr3 {
    qui rdrobust vote margin, vce(`v' state)
    display as text "vce=" %-5s "`v'" "   tau_cl = " %7.4f e(tau_cl) ///
                    "   se_rb = " %7.4f e(se_tau_rb) ///
                    "   pv_rb = " %6.4f e(pv_rb)
}

di _n "Full output for vce(cr3 state):"
rdrobust vote margin, vce(cr3 state)

********************************************************************************
** (2b) Auto-validation: cr* without cluster errors cleanly
********************************************************************************
di _n "Attempt vce(cr3) with no cluster variable (should exit 125):"
cap noi rdrobust vote margin, vce(cr3)
di "rc = " _rc


********************************************************************************
** (3) rdrobustplot -- diagnostic plot
********************************************************************************
di _n "======================================================================"
di       "(3) rdrobustplot: diagnostic plot for a previous rdrobust result"
di       "======================================================================"

qui rdrobust vote margin
rdrobustplot
graph export rdrobustplot_default.png, replace as(png)

qui rdrobust vote margin, vce(cr3 state)
rdrobustplot, shade
graph export rdrobustplot_cr3_shade.png, replace as(png)

qui rdrobust vote margin, covs(termshouse termssenate class)
rdrobustplot, nbins(25 25) title("Senate RD with covariates")
graph export rdrobustplot_covs.png, replace as(png)

di "Saved: rdrobustplot_default.png, rdrobustplot_cr3_shade.png, rdrobustplot_covs.png"

********************************************************************************
** (4) Tables with esttab: combining conventional estimates with RBC inference
********************************************************************************
di _n "======================================================================"
di       "(4) Regression tables via estout/esttab"
di       "======================================================================"
di _n "rdrobust stores three named coefficients in e(b) / e(V):"
di       "    _b[Conventional]    + _se[Conventional]   (p-regression point + SE)"
di       "    _b[Bias-corrected]  + _se[Bias-corrected] (tau_bc with conventional SE;"
di       "                                              kept for display parity, not"
di       "                                              a recommended inference)"
di       "    _b[Robust]          + _se[Robust]         (RBC point + RBC SE)"
di _n "All three rows feed standard Stata tooling (lincom, test, estimates"
di       "table, esttab) without extra scaffolding. The Robust row reproduces"
di       "the author-recommended RBC inference: lincom _b[Robust] returns the"
di       "same CI that rdrobust prints as e(ci_l_rb), e(ci_r_rb) or as the"
di       "third row of the 3x2 matrix e(ci)."

* Install estout from SSC if not already present (one-time per machine).
cap which esttab
if _rc ssc install estout, replace

* Estimate and store four specifications side-by-side.
qui rdrobust vote margin
estimates store m1
qui rdrobust vote margin, covs(class termshouse termssenate)
estimates store m2
qui rdrobust vote margin, vce(cr1 state)
estimates store m3
qui rdrobust vote margin, vce(cr1 state) covs(class termshouse termssenate)
estimates store m4

di _n "-- Table 1: conventional point estimate with RBC confidence interval --"
di       "   (the canonical CCT presentation: report tau_conv, bracket the"
di       "    RBC CI; note that the RBC CI is centered on tau_bc, not tau_conv)"

* Add per-spec scalars/locals used by the LaTeX and console tables below.
* e(ci) is a 3x2 matrix with rows Conventional / Bias-corrected / Robust
* and columns ll / ul; the Robust row is row 3.
foreach m in m1 m2 m3 m4 {
    estimates restore `m'
    local lo : di %9.3f e(ci)[3,1]
    local hi : di %9.3f e(ci)[3,2]
    estadd local RBC_CI = "[`=trim("`lo'")', `=trim("`hi'")']"
    estadd scalar h_bw = e(h_l)
    estadd scalar N_eff = e(N_h_l) + e(N_h_r)
    estimates store `m'
}

esttab m1 m2 m3 m4, ///
    cells(b(fmt(%9.3f))) collabels(none) nostar ///
    keep(Conventional) ///
    scalars(RBC_CI N_eff) ///
    mtitles("sharp" "sharp+covs" "cluster" "cluster+covs") ///
    title("RD estimates: conventional point with RBC CI") ///
    addnotes("Point estimate from p-order local polynomial." ///
             "RBC CI from robust bias correction (centered on tau_bc).")

di _n "-- Table 2: LaTeX export (writes rdrobust_table.tex) --"
esttab m1 m2 m3 m4 using rdrobust_table.tex, replace ///
    cells(b(fmt(%9.3f))) nostar collabels(none) ///
    keep(Conventional) coeflabels(Conventional "RD Effect") ///
    stats(RBC_CI N h_bw N_eff, ///
          labels("Robust 95\% CI" "N" "h" "Eff. N") fmt(0 0 3 0)) ///
    mtitles("sharp" "sharp+covs" "cluster" "cluster+covs") ///
    booktabs nonotes label ///
    title("RD estimates for Senate incumbency")
di "Wrote rdrobust_table.tex"

********************************************************************************
** Done
********************************************************************************
di _n "Done."
