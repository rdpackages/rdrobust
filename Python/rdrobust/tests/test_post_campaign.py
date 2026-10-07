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
    with pytest.raises(ValueError, match="Not enough distinct"):
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


# --- PR #25 review (Matias, 2026-10-01) -----------------------------------

def _review_discrete():
    rng = np.random.default_rng(7)
    x = rng.choice(np.linspace(-1, 1, 12), 60)
    y = 1 + 2 * x + 0.7 * (x >= 0) + rng.normal(size=60)
    return y, x


def test_review2_all_true_rejects_any_undefined_selector():
    # The first row (mserd) is fine here, but msetwo/certwo and their combinations
    # are not identified; the whole table must be checked.
    y, x = _review_discrete()
    with pytest.raises(Exception, match="Not enough variability"):
        quiet(rdbwselect, y, x, all=True)


def test_review2_requested_selector_still_works():
    # Requesting only mserd must not fail because an unrequested selector would.
    y, x = _review_discrete()
    bws = quiet(rdbwselect, y, x).bws.to_numpy().astype(float)
    assert np.all(np.isfinite(bws) & (bws > 0))


def test_review2_all_true_continuous_returns_ten_valid_rows(rd):
    bws = quiet(rdbwselect, rd["y"], rd["x"], all=True).bws
    assert bws.shape[0] == 10
    assert np.all(np.isfinite(bws.to_numpy().astype(float)) & (bws.to_numpy().astype(float) > 0))


def test_review3_few_cluster_warning_ignores_zero_weight_observations():
    i = np.arange(1, 801)
    x = (i - 400.5) / 400
    y = 1 + x + 0.7 * (x >= 0) + 0.2 * np.sin(i * 0.113)
    cluster = i % 40
    weights = (cluster < 4).astype(float)  # 4 contributing clusters per side
    with pytest.warns(UserWarning, match=r"Only 4 \(left\) and 4 \(right\) clusters"):
        with contextlib.redirect_stdout(io.StringIO()):
            rdrobust(y, x, h=0.8, cluster=cluster, weights=weights, vce="cr1")


@pytest.mark.parametrize("name,value", [("h", 0), ("b", [1, 2, 3]), ("h", 2j)])
def test_review4_bandwidth_errors_are_value_errors(name, value):
    x = np.r_[np.arange(-9, 0), np.arange(1, 10)].astype(float)
    y = 1 + 2 * x + 3 * (x >= 0) + 0.2 * np.sin(np.arange(1, 19))
    options = {"h": 12, "b": 15, name: value}
    with pytest.raises(ValueError, match="must contain one or two positive finite numbers"):
        quiet(rdrobust, y, x, **options)


def test_review4_insufficient_support_is_value_error():
    x = np.tile([-2.0, -1.0, 1.0, 2.0], 10)
    y = 1 + 2 * x + 0.2 * np.sin(np.arange(40))
    with pytest.raises(ValueError):
        quiet(rdrobust, y, x, h=25, b=25, vce="hc0")


@pytest.mark.parametrize("rho", [0.123456, 0.00014, 0.00001, 0.5])
def test_review1_b_equals_h_over_rho(rd, rho):
    # Contract b = h/rho for non-round rho (the Stata defect; pinned here too).
    est = quiet(rdrobust, rd["y"], rd["x"], h=0.5, rho=rho)
    assert np.allclose(est.bws.iloc[1].values.astype(float), 0.5 / rho, rtol=1e-12)
    est = quiet(rdrobust, rd["y"], rd["x"], h=[0.4, 0.6], rho=rho)
    assert np.allclose(est.bws.iloc[1].values.astype(float), np.array([0.4, 0.6]) / rho, rtol=1e-12)


# --- independent review of ce435d0 ----------------------------------------

@pytest.mark.parametrize("bwselect", ["msecomb1", "cercomb1"])
def test_combination_selectors_do_not_hide_an_undefined_pilot(bwselect):
    # msesum is undefined here while mserd is fine; min(mserd, msesum) must not
    # return mserd (R and Stata stop on the same data).
    rng = np.random.default_rng(700)
    x = rng.integers(-7, 7, 400).astype(float)
    y = 0.3 * x + (x >= 0) + rng.normal(size=400)
    with pytest.raises(Exception, match="Not enough variability"):
        quiet(rdbwselect, y, x, p=2, bwselect=bwselect)
    bws = quiet(rdbwselect, y, x, p=2).bws.to_numpy().astype(float)
    assert np.all(np.isfinite(bws) & (bws > 0))


def test_p_plus_one_mass_point_clusters_warn_on_degenerate_variance():
    x = np.repeat(np.r_[np.arange(-5, 0), np.arange(0, 5)], 30).astype(float)
    y = 1 + 0.5 * x + (x >= 0.5) + np.random.default_rng(3).normal(size=x.size)
    with pytest.warns(UserWarning, match="estimated variance may be degenerate"):
        with contextlib.redirect_stdout(io.StringIO()):
            est = rdrobust(y, x, c=0.5, h=2, b=5, p=1, cluster=x, vce="cr1")
    assert est.se.iloc[0, 0] < 1e-8


def test_p_plus_one_general_clusters_can_have_positive_variance():
    # Unlike clustering on p+1 mass points, each cluster spans many x values.
    x = np.linspace(-1, 1, 800)
    i = np.arange(x.size)
    cluster = i % 2 + 2 * (x >= 0)
    y = 1 + x + 0.7 * (x >= 0) + 0.3 * (i % 2) + 0.2 * np.sin(i * 0.113)
    with pytest.warns(UserWarning, match="inference may be unreliable") as caught:
        est = rdrobust(y, x, h=0.8, b=0.9, cluster=cluster, vce="cr1")
    assert not any("not identified" in str(w.message) for w in caught)

    # Independent local-linear CR1 calculation. rdrobust uses the union of
    # the h/b windows for the finite-sample correction, with zero h-weights
    # outside the estimation window.
    variances = []
    for side in (x < 0, x >= 0):
        keep = side & (np.abs(x) < 0.9)
        X = np.column_stack([np.ones(keep.sum()), x[keep]])
        w = np.maximum(1 - np.abs(x[keep]) / 0.8, 0)
        invG = np.linalg.inv(X.T @ (w[:, None] * X))
        residual = y[keep] - X @ (invG @ (X.T @ (w * y[keep])))
        scores = np.array([
            X[cluster[keep] == g].T @ (w[cluster[keep] == g] * residual[cluster[keep] == g])
            for g in np.unique(cluster[keep])
        ])
        V = invG @ scores.T @ scores @ invG * ((keep.sum() - 1) / (keep.sum() - 2)) * 2
        variances.append(V[0, 0])
    assert all(v > 0 for v in variances)
    np.testing.assert_allclose(est.se.iloc[0, 0], np.sqrt(sum(variances)), rtol=1e-12)


@pytest.mark.parametrize("bad", [0, 2.5])
def test_nnmatch_is_validated(rd, bad):
    for f in (rdrobust, rdbwselect):
        with pytest.raises(ValueError, match="nnmatch"):
            quiet(f, rd["y"], rd["x"], nnmatch=bad)


def test_rdbwselect_rejects_negative_weights(rd):
    w = np.ones(rd["n"]); w[:5] = -1
    with pytest.raises(ValueError, match="non-negative"):
        quiet(rdbwselect, rd["y"], rd["x"], weights=w)
