# RDROBUST

The package `rdrobust` implements estimation, inference, and graphical procedures for Regression Discontinuity (RD) designs using local polynomial methods.

- `rdrobust`: point estimation and robust bias-corrected inference.
- `rdbwselect`: data-driven bandwidth selection for RD estimation and inference.
- `rdplot`: data-driven RD plots based on binned means and local polynomial fits.

See references for methodological and practical details.

Website: [https://rdpackages.github.io/](https://rdpackages.github.io/).

Source code: [https://github.com/rdpackages/rdrobust](https://github.com/rdpackages/rdrobust).

## Authors
 
Sebastian Calonico (<scalonico@ucdavis.edu>)

Matias D. Cattaneo (<matias.d.cattaneo@gmail.com>)

Max H. Farrell (<mhfarrell@gmail.com>)

Ricardo Masini (<ricardo.masini@gmail.com>)

Rocio Titiunik (<rocio.titiunik@gmail.com>)


## Installation

To install/update use pip:
```
pip install rdrobust
```

## Usage
```python
from rdrobust import rdrobust, rdbwselect, rdplot, plot_rdrobust, rdrobust_RDsenate

# Load bundled Senate dataset (1390 x 14)
df = rdrobust_RDsenate()

# Sharp RD with HC3 heteroskedasticity-robust variance
r = rdrobust(y=df['vote'], x=df['margin'], vce='hc3')

# Cluster-robust variance (CRV3, Pustejovsky-Tipton 2018)
r = rdrobust(y=df['vote'], x=df['margin'], cluster=df['state'], vce='cr3')

# Bandwidth selection
bw = rdbwselect(y=df['vote'], x=df['margin'], all=True)

# Diagnostic plot of an rdrobust object, with effect panel
fig = plot_rdrobust(r, y=df['vote'], x=df['margin'], show_effect=True)
```

- Replication: [rdrobust illustration](https://github.com/rdpackages/rdrobust/blob/main/Python/rdrobust_illustration.py), [rdplot illustration](https://github.com/rdpackages/rdrobust/blob/main/Python/rdplot_illustration.py), [senate data](https://github.com/rdpackages/rdrobust/blob/main/Python/rdrobust_senate.csv).


## Dependencies

- numpy
- pandas
- scipy
- plotnine
- matplotlib

## References

For overviews and introductions, see [rdpackages website](https://rdpackages.github.io).

### Software and Implementation

- Calonico, Cattaneo and Titiunik (2014): [Robust Data-Driven Inference in the Regression-Discontinuity Design](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_Stata.pdf).<br>
_Stata Journal_ 14(4): 909-946.

- Calonico, Cattaneo and Titiunik (2015): [rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_R.pdf).<br>
_R Journal_ 7(1): 38-51.

- Calonico, Cattaneo, Farrell and Titiunik (2017): [rdrobust: Software for Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2017_Stata.pdf).<br>
_Stata Journal_ 17(2): 372-404.

### Technical and Methodological

- Calonico, Cattaneo and Titiunik (2014): [Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA.pdf).<br>
_Econometrica_ 82(6): 2295-2326.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2014_ECMA--Supplemental.pdf).

- Calonico, Cattaneo and Titiunik (2015): [Optimal Data-Driven Regression Discontinuity Plots](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_JASA.pdf).<br>
_Journal of the American Statistical Association_ 110(512): 1753-1769.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Titiunik_2015_JASA--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2018): [On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA.pdf).<br>
_Journal of the American Statistical Association_ 113(522): 767-779.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2018_JASA--Supplement.pdf).

- Calonico, Cattaneo, Farrell and Titiunik (2019): [Regression Discontinuity Designs Using Covariates](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT.pdf).<br>
_Review of Economics and Statistics_ 101(3): 442-451.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell-Titiunik_2019_RESTAT--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2020): [Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ.pdf).<br>
_Econometrics Journal_ 23(2): 192-210.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2020_ECTJ--Supplement.pdf).

- Calonico, Cattaneo and Farrell (2022): [Coverage Error Optimal Confidence Intervals for Local Polynomial Regression](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli.pdf).<br>
_Bernoulli_ 28(4): 2998-3022.<br>
[Supplemental Appendix](https://rdpackages.github.io/references/Calonico-Cattaneo-Farrell_2022_Bernoulli--Supplement.pdf).


<br><br>
