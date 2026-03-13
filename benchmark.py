"""
Benchmark: rdrobust with Numba acceleration vs pure Python.

Tests correctness (numerical equivalence) and measures speedup
at multiple dataset sizes.
"""

import numpy as np
import time
import sys
import os

# Add the source to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), "Python", "rdrobust", "src"))

from rdrobust.funs import (
    _try_numba, rdrobust_res, rdrobust_kweight, nanmat, ncol, qrXXinv, crossprod,
)


def generate_rd_data(n, seed=42):
    """Generate synthetic sharp RD data."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(-1, 1, n)
    y = 3 + 2 * x + 5 * (x >= 0) + rng.normal(0, 1, n)
    return x, y


def generate_rd_data_with_masspoints(n, seed=42):
    """Generate RD data with duplicate running-variable values (mass points)."""
    rng = np.random.default_rng(seed)
    x = np.round(rng.uniform(-1, 1, n), decimals=1)
    y = 3 + 2 * x + 5 * (x >= 0) + rng.normal(0, 1, n)
    return x, y


def prepare_nn_inputs(x, y, T=None, Z=None, c=0):
    """Prepare inputs for rdrobust_res NN path (mimics rdrobust internals)."""
    import pandas as pd

    order = np.argsort(x)
    x = x[order]
    y = y[order]
    T_sorted = T[order] if T is not None else None
    Z_sorted = Z[order] if Z is not None else None

    mask = x < c
    sides = []
    for side_mask in [mask, ~mask]:
        X_s = x[side_mask]
        Y_s = y[side_mask]
        T_s = T_sorted[side_mask] if T_sorted is not None else None
        Z_s = Z_sorted[side_mask] if Z_sorted is not None else None
        n_s = len(X_s)
        aux = pd.DataFrame({"nn": np.ones(n_s), "X": X_s})
        dups = aux.groupby("X")["nn"].transform("sum").values.astype(int)
        dupsid = aux.groupby("X")["nn"].transform("cumsum").values.astype(int)
        sides.append((X_s, Y_s, T_s, Z_s, dups, dupsid, n_s))

    return sides


def rdrobust_res_no_numba(X, y, T, Z, matches, dups, dupsid, d):
    """Force pure Python path by temporarily disabling Numba."""
    import rdrobust.funs as _funs
    saved = _funs._HAS_NUMBA
    _funs._HAS_NUMBA = False
    res = rdrobust_res(X, y, T, Z, 0, 0, "nn", matches, dups, dupsid, d)
    _funs._HAS_NUMBA = saved
    return res


def time_rdrobust_res(X, y, dups, dupsid, n, matches=3, repeats=3):
    """Time rdrobust_res on the NN path."""
    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        res = rdrobust_res(X, y, None, None, 0, 0, "nn", matches, dups, dupsid, 2)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return min(times), res


def time_pure_python_res(X, y, dups, dupsid, n, matches=3, repeats=1):
    """Time the pure Python NN residuals (bypass Numba)."""
    res = nanmat(n, 1)
    t0 = time.perf_counter()
    for pos in range(n):
        rpos = dups[pos] - dupsid[pos]
        lpos = dupsid[pos] - 1
        while lpos + rpos < min(matches, n - 1):
            if pos - lpos - 1 < 0:
                rpos += dups[pos + rpos + 1]
            elif pos + rpos + 1 >= n:
                lpos += dups[pos - lpos - 1]
            elif (X[pos] - X[pos - lpos - 1]) > (X[pos + rpos + 1] - X[pos]):
                rpos += dups[pos + rpos + 1]
            elif (X[pos] - X[pos - lpos - 1]) < (X[pos + rpos + 1] - X[pos]):
                lpos += dups[pos - lpos - 1]
            else:
                rpos += dups[pos + rpos + 1]
                lpos += dups[pos - lpos - 1]
        ind_J = np.arange(max(0, pos - lpos), min(n, pos + rpos) + 1)
        y_J = sum(y[ind_J]) - y[pos]
        Ji = len(ind_J) - 1
        res[pos, 0] = np.sqrt(Ji / (Ji + 1)) * (y[pos] - y_J / Ji)
    t1 = time.perf_counter()
    return t1 - t0, res


def run_full_rdrobust(x, y, repeats=3):
    """Time a full rdrobust() call end-to-end."""
    from rdrobust.rdrobust import rdrobust as rd

    times = []
    for _ in range(repeats):
        t0 = time.perf_counter()
        result = rd(y, x)
        t1 = time.perf_counter()
        times.append(t1 - t0)
    return min(times), result


def check_correctness(label, res_python, res_numba):
    """Compare Python and Numba results, print pass/fail."""
    max_diff = np.max(np.abs(res_python - res_numba))
    print(f"   Max absolute difference: {max_diff:.2e}")
    if max_diff < 1e-10:
        print("   PASS\n")
        return True
    print("   FAIL\n")
    return False


def main():
    has_numba = _try_numba()
    print(f"Numba available: {has_numba}")

    if has_numba:
        print("Numba detected - JIT compilation on first call\n")

        # ── Warm up Numba JIT ──
        print("Warming up Numba JIT...")
        x_warm, y_warm = generate_rd_data(1000, seed=0)
        sides = prepare_nn_inputs(x_warm, y_warm)
        X_s, Y_s, _, _, dups, dupsid, n_s = sides[0]
        _ = rdrobust_res(X_s, Y_s, None, None, 0, 0, "nn", 3, dups, dupsid, 2)
        print("JIT warm-up done.\n")

    # ── Correctness checks ──
    print("=" * 60)
    print("CORRECTNESS CHECKS")
    print("=" * 60)

    # 1. Sharp RD, continuous X
    print("\n1. Sharp RD, continuous X (n=5,000)")
    x_check, y_check = generate_rd_data(5000, seed=1)
    sides = prepare_nn_inputs(x_check, y_check)
    X_s, Y_s, _, _, dups, dupsid, n_s = sides[0]
    _, res_python = time_pure_python_res(X_s, Y_s, dups, dupsid, n_s)
    _, res_numba = time_rdrobust_res(X_s, Y_s, dups, dupsid, n_s)
    if not check_correctness("sharp", res_python, res_numba):
        return

    # 2. Sharp RD with mass points (duplicate X values)
    print("2. Sharp RD, mass points (n=5,000)")
    x_mp, y_mp = generate_rd_data_with_masspoints(5000, seed=1)
    sides_mp = prepare_nn_inputs(x_mp, y_mp)
    X_s, Y_s, _, _, dups, dupsid, n_s = sides_mp[0]
    _, res_py_mp = time_pure_python_res(X_s, Y_s, dups, dupsid, n_s)
    _, res_nb_mp = time_rdrobust_res(X_s, Y_s, dups, dupsid, n_s)
    if not check_correctness("mass points", res_py_mp, res_nb_mp):
        return

    # 3. Fuzzy RD (y + T)
    print("3. Fuzzy RD (n=5,000)")
    rng = np.random.default_rng(1)
    x_fz = rng.uniform(-1, 1, 5000)
    T_fz = ((x_fz >= 0).astype(float) + rng.binomial(1, 0.1, 5000)).clip(0, 1)
    y_fz = 3 + 2 * x_fz + 1.5 * T_fz + rng.normal(0, 1, 5000)
    sides_fz = prepare_nn_inputs(x_fz, y_fz, T=T_fz)
    X_s, Y_s, T_s, _, dups, dupsid, n_s = sides_fz[0]
    res_py_fz = rdrobust_res_no_numba(X_s, Y_s, T_s, None, 3, dups, dupsid, 2)
    res_nb_fz = rdrobust_res(X_s, Y_s, T_s, None, 0, 0, "nn", 3, dups, dupsid, 2)
    if not check_correctness("fuzzy", res_py_fz, res_nb_fz):
        return

    # 4. With covariates (y + T + Z)
    print("4. With covariates (n=5,000)")
    Z_fz = rng.normal(0, 1, (5000, 2))
    sides_cov = prepare_nn_inputs(x_fz, y_fz, T=T_fz, Z=Z_fz)
    X_s, Y_s, T_s, Z_s, dups, dupsid, n_s = sides_cov[0]
    res_py_cov = rdrobust_res_no_numba(X_s, Y_s, T_s, Z_s, 3, dups, dupsid, 2)
    res_nb_cov = rdrobust_res(X_s, Y_s, T_s, Z_s, 0, 0, "nn", 3, dups, dupsid, 2)
    if not check_correctness("covariates", res_py_cov, res_nb_cov):
        return

    # ── Benchmark rdrobust_res NN path ──
    print("=" * 60)
    print("BENCHMARK: rdrobust_res (NN path, left side only)")
    print("=" * 60)
    print(f"{'n':>10} {'Python (s)':>12} {'Numba (s)':>12} {'Speedup':>10}")
    print("-" * 50)

    for n in [5_000, 20_000, 50_000, 100_000, 200_000, 500_000]:
        x_bench, y_bench = generate_rd_data(n * 2, seed=42)
        sides = prepare_nn_inputs(x_bench, y_bench)
        X_s, Y_s, _, _, dups, dupsid, n_s = sides[0]

        # Numba timing
        t_numba, _ = time_rdrobust_res(X_s, Y_s, dups, dupsid, n_s, repeats=3)

        # Python timing (skip for large n - too slow)
        if n <= 50_000:
            t_python, _ = time_pure_python_res(X_s, Y_s, dups, dupsid, n_s, repeats=1)
            speedup = t_python / t_numba
            print(f"{n_s:>10,} {t_python:>12.4f} {t_numba:>12.4f} {speedup:>9.0f}x")
        else:
            # Extrapolate Python time (linear in n)
            print(f"{n_s:>10,} {'(skipped)':>12} {t_numba:>12.4f} {'~':>9}")

    # ── End-to-end benchmark ──
    print("\n" + "=" * 60)
    print("BENCHMARK: Full rdrobust() end-to-end")
    print("=" * 60)
    print(f"{'n':>10} {'Time (s)':>12}")
    print("-" * 25)

    for n in [10_000, 50_000, 100_000, 500_000]:
        x_e2e, y_e2e = generate_rd_data(n, seed=42)
        t_e2e, result = run_full_rdrobust(x_e2e, y_e2e, repeats=2)
        coef = result.coef.iloc[0, 0]
        print(f"{n:>10,} {t_e2e:>12.3f}   tau={coef:.4f}")


if __name__ == "__main__":
    main()
