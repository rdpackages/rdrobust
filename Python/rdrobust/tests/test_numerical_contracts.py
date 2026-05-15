import contextlib
import io

import numpy as np

from rdrobust import rdbwselect, rdrobust
from rdrobust.datasets import rdrobust_RDsenate


def senate_covariate_data():
    data = rdrobust_RDsenate()
    return (
        data["vote"],
        data["margin"],
        data[["termshouse", "termssenate", "population"]],
    )


def test_covariate_bandwidths_match_type2_iqr_baseline():
    y, x, covs = senate_covariate_data()
    result = rdbwselect(y, x, covs=covs, vce="hc1")

    np.testing.assert_allclose(
        result.bws.loc["mserd"].to_numpy(dtype=float),
        np.array([17.976890, 17.976890, 28.963712, 28.963712]),
        rtol=1e-7,
        atol=1e-7,
    )


def test_manual_bandwidth_masspoint_counts_match_r_baseline():
    y, x, covs = senate_covariate_data()
    result = rdrobust(y, x, covs=covs, h=[10, 10], b=[20, 20], vce="hc1")

    assert result.N == [491, 617]
    assert result.M == [491, 580]
    assert list(result.N_h) == [215, 181]
    assert list(result.N_b) == [336, 300]


def test_repr_returns_text_without_printing():
    y, x, _ = senate_covariate_data()
    result = rdrobust(y, x)

    stream = io.StringIO()
    with contextlib.redirect_stdout(stream):
        text = repr(result)

    assert stream.getvalue() == ""
    assert "Call: rdrobust" in text
    assert "Sharp RD estimates" in text
