"""Defects found after the 2026-08 fix campaign (D1-D3, R-17), mirrored from R."""

import io
import contextlib
import warnings

import numpy as np
import pytest

from rdrobust import rdbwselect, rdrobust


def quiet(f, *a, **k):
    with contextlib.redirect_stdout(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return f(*a, **k)


@pytest.fixture
def rd():
    rng = np.random.default_rng(7)
    n = 1500
    x = rng.uniform(-1, 1, n)
    y = 1 + 2 * x + 0.7 * (x >= 0) + rng.normal(0, 1, n)
    return {"y": y, "x": x, "n": n}


def test_d2_near_collinear_covariate_on_one_side(rd):
    # z2 equals z1 up to 1e-9 noise on the left only. The selected bandwidth
    # must not depend on that noise.
    rng = np.random.default_rng(1)
    x, n = rd["x"], rd["n"]
    z1, zr = rng.normal(size=n), rng.normal(size=n)
    y = 1 + x + 0.5 * (x >= 0) + 0.3 * z1 + rng.normal(size=n)
    hs = []
    for s in range(4):
        e = np.random.default_rng(100 + s).normal(size=n)
        z2 = np.where(x < 0, z1 + 1e-9 * e, zr)
        hs.append(quiet(rdbwselect, y, x, covs=np.column_stack([z1, z2])).bws.iloc[0, 0])
    assert max(hs) - min(hs) < 1e-6


def test_d3_few_clusters_warn(rd):
    few = np.digitize(rd["x"], np.linspace(-1, 1, 9)[1:-1])  # 4 per side
    with pytest.warns(UserWarning, match="clusters"):
        with contextlib.redirect_stdout(io.StringIO()):
            rdrobust(rd["y"], rd["x"], h=0.5, cluster=few, vce="cr1")


def test_d3_many_clusters_do_not_warn(rd):
    many = np.arange(rd["n"]) % 40
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        with contextlib.redirect_stdout(io.StringIO()):
            rdrobust(rd["y"], rd["x"], h=0.5, cluster=many, vce="cr1")


def test_r17_single_observation_side_fails_informatively(rd):
    rng = np.random.default_rng(11)
    x = rd["x"]
    d = (rng.uniform(size=rd["n"]) < np.where(x >= 0, 0.85, 0.2)).astype(float)
    keep = (x >= 0) | (np.arange(rd["n"]) == np.flatnonzero(x < 0)[0])
    with pytest.raises(Exception, match="Not enough distinct"):
        quiet(rdrobust, rd["y"][keep], x[keep], fuzzy=d[keep])


@pytest.mark.parametrize("K", range(3, 10))
@pytest.mark.parametrize("p", [1, 2])
def test_discrete_running_variable_gives_bandwidth_or_informative_error(K, p):
    for seed in (1, 2):
        rng = np.random.default_rng(seed)
        x = rng.choice(np.linspace(-1, 1, K), 2000)
        y = 1 + 2 * x + 0.7 * (x >= 0) + rng.normal(size=2000)
        try:
            bws = quiet(rdbwselect, y, x, p=p).bws.iloc[0].values.astype(float)
        except Exception as e:
            assert "Not enough" in str(e), f"K={K} seed={seed} p={p}: {e!r}"
        else:
            assert np.all(np.isfinite(bws) & (bws > 0)), f"K={K} seed={seed} p={p}"
