# -*- coding: utf-8 -*-
"""Bundled datasets for rdrobust."""
from importlib import resources
import pandas as pd


def rdrobust_RDsenate():
    """Load the RDsenate dataset.

    A data frame with 1390 observations on 17 variables from Cattaneo,
    Frandsen, and Titiunik (2015), "Randomization Inference in the
    Regression Discontinuity Design: An Application to the Study of Party
    Advantages in the U.S. Senate," Journal of Causal Inference 3(1): 1-24.

    Columns
    -------
    state : str
        U.S. state name; usable as a cluster variable.
    year : int
        Election year.
    margin : float
        Democratic vote margin at election t (running variable).
    vote : float
        Democratic vote share at election t+2 (main outcome).
    class_, termshouse, termssenate : int
        Senate class and number of terms served in House/Senate.
    dopen, dmidterm, dpresdem, demwinprv1, demwinprv2 : int
        Binary indicators.
    population, presdemvoteshlag1, demvoteshlag1, demvoteshlag2,
    demvoteshfor1 : float
        Additional numeric covariates.
    """
    with resources.files(__package__).joinpath("data/rdrobust_RDsenate.csv").open("r") as f:
        return pd.read_csv(f)
