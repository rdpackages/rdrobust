********************************************************************************
** RDROBUST Stata Package
** Do-file for RDPLOT Illustration
**
** Shows how to construct RD plots manually from the output of rdplot.
** rdplot with the `genvars' option writes the bin means, pointwise CIs, and
** per-observation fitted values back to the dataset; `hide' suppresses the
** default plot. With these variables you can build any plot you want.
** Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell and Rocio Titiunik
********************************************************************************
** net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/master/stata) replace
********************************************************************************
clear all

********************************************************************************
** Load data and run rdplot with `genvars, hide, ci()' so we get the bins and
** pointwise CIs back in the data without the default plot.
********************************************************************************
findfile rdrobust_senate.dta
use "`r(fn)'", clear
global y vote
global x margin
global c 0
su $x
global x_min = r(min)
global x_max = r(max)

rdplot $y $x, genvars hide ci(95)

** Polynomial order actually used by rdplot (defaults to 4, but honors p()).
local p = e(p)

** Shared plot elements (reused below):
**   `scatter'    -- bin means dot layer
**   `fit_expr'   -- left/right polynomial fits via `twoway function' on e(eq_*)
**   `fit_hat'    -- left/right polynomial fits via rdplot_hat_y (no eval needed)
**   `common'    -- axis, legend, and title options
local scatter  (scatter rdplot_mean_y rdplot_mean_bin, sort msize(small) mcolor(gs10))
local fit_expr (function `e(eq_l)', range($x_min $c) lcolor(black) sort lwidth(medthin) lpattern(solid)) ///
               (function `e(eq_r)', range($c $x_max) lcolor(black) sort lwidth(medthin) lpattern(solid))
local fit_hat  (line rdplot_hat_y $x if $x< $c, sort lcolor(black) lwidth(medthin)) ///
               (line rdplot_hat_y $x if $x>=$c, sort lcolor(black) lwidth(medthin))
local common xline($c, lcolor(black) lwidth(medthin)) xscale(r($x_min $x_max)) ///
             xtitle("Running variable") ytitle("Outcome") ///
             title("Regression function fit", color(gs0))

********************************************************************************
** (1) Default RDPLOT: bin means + polynomial fit via e(eq_l)/e(eq_r).
********************************************************************************
twoway `scatter' `fit_expr', `common' ///
    legend(pos(6) cols(2) order(1 "Sample average within bin" ///
                                2 "Polynomial fit of order `p'"))

********************************************************************************
** (2) Same plot built with rdplot_hat_y instead of e(eq_*).
**     Uses the per-observation fitted values stored by `genvars'; no
**     polynomial-string eval and no dependency on the Mata coefficient matrices.
********************************************************************************
twoway `scatter' `fit_hat', `common' ///
    legend(pos(6) cols(2) order(1 "Sample average within bin" ///
                                2 "Polynomial fit of order `p'"))

********************************************************************************
** (3) RDPLOT with pointwise 95% confidence intervals as caps.
********************************************************************************
twoway (rcap rdplot_ci_l rdplot_ci_r rdplot_mean_bin, color(gs11)) ///
       `scatter' `fit_expr', `common' ///
    legend(pos(6) cols(3) order(2 "Sample average within bin" ///
                                3 "Polynomial fit of order `p'" ///
                                1 "95% CI"))

********************************************************************************
** (4) RDPLOT with shaded pointwise 95% CI bands (one rarea per side).
********************************************************************************
twoway (rarea rdplot_ci_l rdplot_ci_r rdplot_mean_bin if rdplot_id<0, sort color(gs11)) ///
       (rarea rdplot_ci_l rdplot_ci_r rdplot_mean_bin if rdplot_id>0, sort color(gs11)) ///
       `scatter' `fit_expr', `common' ///
    legend(pos(6) cols(3) order(3 "Sample average within bin" ///
                                4 "Polynomial fit of order `p'" ///
                                1 "95% CI"))
