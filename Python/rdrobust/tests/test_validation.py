"""Input validation, ported from R's precheck block (PY-11).

Every case here used to die inside the estimator with a message naming linear
algebra or an unbound local rather than the argument at fault, or -- worse --
named the wrong argument entirely. R rejects all of them up front; these tests
pin the parity.
"""

import numpy as np
import pytest

from rdrobust import rdbwselect, rdrobust


@pytest.fixture
def rd():
    rng = np.random.default_rng(3)
    n = 800
    x = rng.uniform(-1, 1, n)
    y = 1 + 2 * x + 0.7 * (x >= 0) + rng.normal(0, 1, n)
    return {"y": y, "x": x, "n": n}


# --- bandwidths -----------------------------------------------------------
# A bad h or b used to reach the Cholesky factorization as a degenerate
# design ("LinAlgError: Internal potrf return info = [1]"), and a length-3 h
# matched neither the scalar nor the length-2 branch, so the estimator ran on
# with h_l never assigned ("UnboundLocalError: cannot access local variable
# 'h_l'") -- an internal-logic error surfacing as the user's error message.

@pytest.mark.parametrize("bad_h", [-0.5, 0.0, np.nan, np.inf, [0.3, 0.4, 0.5]])
def test_rdrobust_rejects_bad_h(rd, bad_h):
    with pytest.raises(Exception, match="h must be a positive scalar"):
        rdrobust(rd["y"], rd["x"], h=bad_h)


@pytest.mark.parametrize("bad_b", [-0.8, 0.0, np.nan, [0.3, 0.4, 0.5]])
def test_rdrobust_rejects_bad_b(rd, bad_b):
    with pytest.raises(Exception, match="b must be a positive scalar"):
        rdrobust(rd["y"], rd["x"], h=0.5, b=bad_b)


def test_valid_bandwidth_shapes_still_accepted(rd):
    assert rdrobust(rd["y"], rd["x"], h=0.4).coef is not None
    assert rdrobust(rd["y"], rd["x"], h=[0.4, 0.5]).coef is not None
    assert rdrobust(rd["y"], rd["x"], h=0.4, b=0.7).coef is not None


# --- weights --------------------------------------------------------------

def test_all_zero_weights_are_rejected(rd):
    with pytest.raises(Exception, match="at least one positive value"):
        rdrobust(rd["y"], rd["x"], weights=np.zeros(rd["n"]))


def test_negative_weights_are_named_as_the_problem(rd):
    # Negative weights were swept into the NA mask and dropped silently. With
    # the negatives all on one side that emptied that side of the cutoff, and
    # the reported error was "c should be set within the range of x" -- an
    # error about the cutoff, for a mistake in weights.
    w = np.where(rd["x"] < 0, -1.0, 1.0)
    with pytest.raises(Exception, match="must be non-negative"):
        rdrobust(rd["y"], rd["x"], weights=w)


# --- cutoff outside the support ------------------------------------------

def test_rdrobust_rejects_cutoff_outside_support(rd):
    with pytest.raises(Exception, match="within the range of x"):
        rdrobust(rd["y"], rd["x"], c=5.0)


def test_rdbwselect_rejects_cutoff_outside_support(rd):
    # rdbwselect had no check at all: one side came out empty and the failure
    # was "zero-size array to reduction operation minimum".
    with pytest.raises(Exception, match="within the range of x"):
        rdbwselect(rd["y"], rd["x"], c=5.0)


def test_cutoff_on_the_boundary_is_rejected(rd):
    with pytest.raises(Exception, match="within the range of x"):
        rdrobust(rd["y"], rd["x"], c=float(np.max(rd["x"])))
