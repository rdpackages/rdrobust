# Robust Local Polynomial Methods for RD designs

The package `rdrobust` implements estimation, inference, and graphical procedures for Regression Discontinuity (RD) designs using local polynomial methods.

- `rdrobust`: point estimation and robust bias-corrected inference.
- `rdbwselect`: data-driven bandwidth selection for RD estimation and inference.
- `rdplot`: data-driven RD plots based on binned means and local polynomial fits.
- `rdrobustplot`: postestimation diagnostic plot for Stata.


## Python Implementation

To install/update in Python type:
```
pip install rdrobust
```

- Help: [PYPI repository](https://pypi.org/project/rdrobust/).

- Replication: [rdrobust illustration](Python/rdrobust_illustration.py), [rdplot illustration](Python/rdplot_illustration.py), [senate data](Python/rdrobust_senate.csv).

## R Implementation

To install/update in R type:
```
install.packages('rdrobust')
```

- Help: [R Manual](https://cran.r-project.org/web/packages/rdrobust/rdrobust.pdf), [CRAN repository](https://cran.r-project.org/package=rdrobust).

- Examples/data: [rdrobust illustration](R/rdrobust_illustration.r), [rdplot illustration](R/rdplot_illustration.R), [senate data](R/rdrobust_senate.csv).

## Stata Implementation

To install/update in Stata type:
```
net install rdrobust, from(https://raw.githubusercontent.com/rdpackages/rdrobust/main/stata) replace
```

- Help: [rdrobust](stata/rdrobust.pdf), [rdbwselect](stata/rdbwselect.pdf), [rdplot](stata/rdplot.pdf), [rdrobustplot](stata/rdrobustplot.pdf).

- Replication: [rdrobust illustration](stata/rdrobust_illustration.do), [rdplot illustration](stata/rdplot_illustration.do), [new features illustration](stata/rdrobust_illustration_new.do), [senate data](stata/rdrobust_senate.dta).


## References

For overviews and introductions, see the [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Calonico, Cattaneo and Titiunik (2014): [Robust Data-Driven Inference in the Regression-Discontinuity Design](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf). _Stata Journal_ 14(4): 909-946.
- Calonico, Cattaneo and Titiunik (2015): [rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_R.pdf). _R Journal_ 7(1): 38-51.
- Calonico, Cattaneo, Farrell and Titiunik (2017): [rdrobust: Software for Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf). _Stata Journal_ 17(2): 372-404.

### Technical and Methodological

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf). _Econometrica_ 82(6): 2295-2326. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).
- Calonico, Cattaneo and Titiunik (2015): [Optimal Data-Driven Regression Discontinuity Plots](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_JASA.pdf). _Journal of the American Statistical Association_ 110(512): 1753-1769. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_JASA--Supplement.pdf).
- Calonico, Cattaneo and Farrell (2018): [On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf). _Journal of the American Statistical Association_ 113(522): 767-779. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA--Supplement.pdf).
- Calonico, Cattaneo, Farrell and Titiunik (2019): [Regression Discontinuity Designs Using Covariates](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf). _Review of Economics and Statistics_ 101(3): 442-451. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf).
- Calonico, Cattaneo and Farrell (2020): [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf). _Econometrics Journal_ 23(2): 192-210. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf).
- Calonico, Cattaneo and Farrell (2022): [Coverage Error Optimal Confidence Intervals for Local Polynomial Regression](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf). _Bernoulli_ 28(4): 2998-3022. [Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli--Supplement.pdf).

## Funding

This work was supported in part by the National Science Foundation through grants [SES-1357561](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1357561), [SES-1459931](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1459931), [SES-1947805](https://www.nsf.gov/awardsearch/showAward?AWD_ID=1947805), and [SES-2019432](https://www.nsf.gov/awardsearch/showAward?AWD_ID=2019432).
