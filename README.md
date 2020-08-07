# RDROBUST

The **rdrobust** package provides <font face="courier new, monospace">Stata</font> and <font face="courier new, monospace">R</font> implementations of statistical inference and graphical procedures for Regression Discontinuity designs employing local polynomial and partitioning methods. It provides point estimators, confidence intervals estimators, bandwidth selectors, automatic RD plots, and other related features.

This work was supported in part by the National Science Foundation through grants <a href="http://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561" target="_blank">SES-1357561</a>, <a href="http://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931" target="_blank">SES-1459931</a>, and <a href="http://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805" target="_blank">SES-1947805</a>.

<br>

## Major Upgrades

This package was first released in Spring 2014, and had two major upgrades in Fall 2016 and in Winter 2020.

- _Fall 2016 new features include_: (i) major speed improvements; (ii) covariate-adjusted bandwidth selection, point estimation, and robust inference; (iii) cluster-robust bandwidth selection, point estimation, and robust inference; (iv) weighted global polynomial fits and pointwise confidence bands for RD plots; and (v) several new bandwidths selectors (e.g., different bandwidths for control and treatment groups, coverage error optimal bandwidths, and optimal bandwidths for fuzzy designs).

- _Winter 2020 new features include_: (i) discrete running variable checks and adjustments; (ii) bandwidth selection adjustments for too few mass points in and/or overshooting of the support of the running variable; (iii) RD Plots with additional covariates plotted at their mean (previously the package set additional covariates at zero); (iv) automatic removal of co-linear additional covariates; (v) turn on/off standardization of variables (which may lead to small numerical/rounding discrepancies with prior versions); and (vi) <font face="courier new, monospace">rdplot</font> output using <font face="courier new, monospace">ggplot2</font> in <font face="courier new, monospace">R</font>.

<br>

## References

### Technical and Methodological:

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf), _Econometrica_ 82(6): 2295-2326. [Supplemental Appendix](references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).

- Calonico, Cattaneo and Titiunik (2015): (Optimal Data-Driven Regression Discontinuity Plots)[references/Calonico-Cattaneo-Titiunik_2015_JASA.pdf], _Journal of the American Statistical Association_ 110(512): 1753-1769. [Supplemental Appendix](references/Calonico-Cattaneo-Titiunik_2015_JASA--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2018): <a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell_2018_JASA.pdf?attredirects=0" target="_blank">On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference</a>, <i>Journal of the American Statistical Association</i> 113(522): 767-779. [<a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell_2018_JASA--Supplement.pdf?attredirects=0" target="_blank">Supplemental Appendix</a>]</li>

- Calonico, Cattaneo, Farrell and Titiunik (2019): <a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf?attredirects=0" target="_blank">Regression Discontinuity Designs Using Covariates</a>, <i>Review of Economics and Statistics</i> 101(3): 442-451. [<a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf?attredirects=0" target="_blank">Supplemental Appendix</a>]

- Calonico, Cattaneo and Farrell (2020): <a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf?attredirects=0" target="_blank">Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs</a>, <i>Econometrics Journal</i> 23(2): 192-210. [<a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf?attredirects=0" target="_blank">Supplemental Appendix</a>]

- Calonico, Cattaneo and Farrell (2020): <a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell_2020_CEopt.pdf?attredirects=0" target="_blank">Coverage Error Optimal Confidence Intervals for Local Polynomial Regression</a>. [<a href="https://sites.google.com/site/rdpackages/rdrobust/Calonico-Cattaneo-Farrell_2020_CEopt--Supplement.pdf?attredirects=0" target="_blank">Supplemental Appendix</a>]

### Software and Implementation

- Calonico, Cattaneo and Titiunik (2014): [Robust Data-Driven Inference in the Regression-Discontinuity Design](references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf), _Stata Journal_ 14(4): 909-946.

- Calonico, Cattaneo and Titiunik (2015): [rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs](references/Calonico-Cattaneo-Titiunik_2015_R.pdf), _R Journal_ 7(1): 38-51.

- Calonico, Cattaneo, Farrell and Titiunik (2017): [rdrobust: Software for Regression Discontinuity Designs](references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf), _Stata Journal_ 17(2): 372-404.

<br>

## Stata Implementation:</b><br>

To install/update in Stata type:
```
net install rdrobust, from(https://github.com/rdpackages/rdrobust/stata) replace
```

<ul>
<li style="margin-bottom:10px">Help files: <a href="https://sites.google.com/site/rdpackages/rdrobust/stata/rdrobust.pdf?attredirects=0" target="_blank"><font face="courier new, monospace">rdrobust</font></a>, <a href="https://sites.google.com/site/rdpackages/rdrobust/stata/rdbwselect.pdf?attredirects=0" target="_blank"><font face="courier new, monospace">rdbwselect</font></a>, <font face="courier new, monospace"><a href="https://sites.google.com/site/rdpackages/rdrobust/stata/rdplot.pdf?attredirects=0" target="_blank">rdplot</a></font> -- Replication files: <a href="https://sites.google.com/site/rdpackages/rdrobust/stata/rdrobust_illustration.do?attredirects=0" target="_blank">do-file</a>, <a href="https://sites.google.com/site/rdpackages/rdrobust/stata/rdplot_illustration.do?attredirects=0" target="_blank"><font face="courier new, monospace">rdplot</font> illustration</a>, <a href="https://sites.google.com/site/rdpackages/rdrobust/stata/rdrobust_senate.dta?attredirects=0" target="_blank">senate data</a></li>

<li style="margin-bottom:10px">Repository for manual installation: <a href="https://sites.google.com/site/rdpackages/rdrobust/stata" target="_blank">https://sites.google.com/site/rdpackages/rdrobust/stata</a></li>
</ul>
<br>

<b>Implementation in R:</b><br>

<ul><li>To install/update in R type:</li></ul>
<font face="courier new, monospace" style="margin-left:60px">install.packages('rdrobust')</font>

<ul>
<li style="margin-bottom:10px"><a href="https://cran.r-project.org/web/packages/rdrobust/rdrobust.pdf" target="_blank">Manual</a> -- Replication files: <a href="https://sites.google.com/site/rdpackages/rdrobust/r/rdrobust_illustration.r?attredirects=0" target="_blank">R-script</a>, <a href="https://sites.google.com/site/rdpackages/rdrobust/r/rdplot_illustration.R?attredirects=0" target="_blank"><font face="courier new, monospace">rdplot</font> illustration</a>, <a href="https://sites.google.com/site/rdpackages/rdrobust/r/rdrobust_senate.csv?attredirects=0" target="_blank">senate data</a></li>
<li style="margin-bottom:10px"><a href="https://cran.r-project.org/package=rdrobust" target="_blank">CRAN repository</a> -- Repository for manual installation: <a href="https://sites.google.com/site/rdpackages/rdrobust/r" target="_blank">https://sites.google.com/site/rdpackages/rdrobust/r</a></li>
</ul>
<br>

<b>Last update</b>: June 4, 2020.

<br>
<br>
<br>
