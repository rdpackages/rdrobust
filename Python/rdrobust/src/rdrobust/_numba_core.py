"""
Numba-accelerated nearest-neighbor residual computation for rdrobust.

Replaces the O(n) Python loop in the NN VCE path with JIT-compiled code.
The NN residual computation is the dominant bottleneck in rdrobust, called
7-19 times per estimation. Each call loops over all n observations in
pure Python. This module provides ~100x speedup on that hot path.
"""

import numpy as np
from numba import njit, prange


@njit(cache=True)
def nn_residuals(X, D, matches, dups, dupsid, n, ncols):
    """
    Compute nearest-neighbor residuals for all data columns.

    Parameters
    ----------
    X : float64 array (n,)
        Running variable, sorted.
    D : float64 array (n, ncols)
        Stacked data columns: [y] or [y, T] or [y, T, Z1, Z2, ...].
    matches : int
        Minimum number of neighbors (nnmatch, default 3).
    dups : int64 array (n,)
        Number of duplicate X values at each position.
    dupsid : int64 array (n,)
        Cumulative duplicate ID at each position.
    n : int
        Number of observations.
    ncols : int
        Number of columns in D.

    Returns
    -------
    res : float64 array (n, ncols)
    """
    res = np.empty((n, ncols))
    target = min(matches, n - 1)

    for pos in range(n):
        rpos = dups[pos] - dupsid[pos]
        lpos = dupsid[pos] - 1

        while lpos + rpos < target:
            left_ok = pos - lpos - 1 >= 0
            right_ok = pos + rpos + 1 < n

            if not left_ok:
                rpos += dups[pos + rpos + 1]
            elif not right_ok:
                lpos += dups[pos - lpos - 1]
            else:
                d_left = X[pos] - X[pos - lpos - 1]
                d_right = X[pos + rpos + 1] - X[pos]
                if d_left > d_right:
                    rpos += dups[pos + rpos + 1]
                elif d_left < d_right:
                    lpos += dups[pos - lpos - 1]
                else:
                    rpos += dups[pos + rpos + 1]
                    lpos += dups[pos - lpos - 1]

        lo = max(0, pos - lpos)
        hi = min(n, pos + rpos) + 1  # exclusive
        Ji = (hi - lo) - 1
        scale = np.sqrt(Ji / (Ji + 1.0))

        for col in range(ncols):
            col_sum = 0.0
            for k in range(lo, hi):
                col_sum += D[k, col]
            col_sum -= D[pos, col]
            res[pos, col] = scale * (D[pos, col] - col_sum / Ji)

    return res


@njit(cache=True, parallel=True)
def nn_residuals_parallel(X, D, matches, dups, dupsid, n, ncols):
    """Parallel version of nn_residuals using prange."""
    res = np.empty((n, ncols))
    target = min(matches, n - 1)

    for pos in prange(n):
        rpos = dups[pos] - dupsid[pos]
        lpos = dupsid[pos] - 1

        while lpos + rpos < target:
            left_ok = pos - lpos - 1 >= 0
            right_ok = pos + rpos + 1 < n

            if not left_ok:
                rpos += dups[pos + rpos + 1]
            elif not right_ok:
                lpos += dups[pos - lpos - 1]
            else:
                d_left = X[pos] - X[pos - lpos - 1]
                d_right = X[pos + rpos + 1] - X[pos]
                if d_left > d_right:
                    rpos += dups[pos + rpos + 1]
                elif d_left < d_right:
                    lpos += dups[pos - lpos - 1]
                else:
                    rpos += dups[pos + rpos + 1]
                    lpos += dups[pos - lpos - 1]

        lo = max(0, pos - lpos)
        hi = min(n, pos + rpos) + 1
        Ji = (hi - lo) - 1
        scale = np.sqrt(Ji / (Ji + 1.0))

        for col in range(ncols):
            col_sum = 0.0
            for k in range(lo, hi):
                col_sum += D[k, col]
            col_sum -= D[pos, col]
            res[pos, col] = scale * (D[pos, col] - col_sum / Ji)

    return res


# Threshold for switching to parallel version
PARALLEL_THRESHOLD = 50_000
