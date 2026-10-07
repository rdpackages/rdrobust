import contextlib
import io

import numpy as np

from rdrobust import rdbwselect, rdplot, rdrobust
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
def test_bin_edges_stay_attached_to_their_own_bin_when_a_side_has_holes():
    # The left support has a gap in [-0.55, -0.25], so several evenly-spaced
    # left bins come out empty while the right side is fully occupied. Empty
    # bins are dropped from vars_bins, and the surviving rows must keep the
    # edges of the bins they were actually computed from.
    rng = np.random.default_rng(20260813)
    x = np.concatenate([
        rng.uniform(-1, -0.55, 200),
        rng.uniform(-0.25, 0, 200),
        rng.uniform(0, 1, 400),
    ])
    y = 3 + 2 * x + 4 * (x >= 0) + rng.normal(0, 0.3, len(x))

    out = rdplot(y, x, nbins=[12, 12], binselect="es", hide=True)
    vb = out.vars_bins

    # The gap must actually have emptied some left bins, or the test is vacuous.
    assert int((vb.rdplot_mean_x < 0).sum()) < out.J[0]

    # Each bin's mean falls inside that bin, and the reported edges bracket the
    # midpoint that rdplot computed independently.
    assert (vb.rdplot_mean_x >= vb.rdplot_min_bin).all()
    assert (vb.rdplot_mean_x <= vb.rdplot_max_bin).all()
    np.testing.assert_allclose(
        (vb.rdplot_min_bin + vb.rdplot_max_bin) / 2,
        vb.rdplot_mean_bin,
        rtol=0, atol=1e-12,
    )
