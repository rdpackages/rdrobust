#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Rdrobust module

Created on Sat Jun  5 14:07:58 2021

@author: rmasini
"""

import numpy as np
import pandas as pd
import math
import io
from contextlib import redirect_stdout
from scipy.linalg import qr, cho_factor, cho_solve

try:
    from numba import njit
    _NB_AVAILABLE = True
except ImportError:
    _NB_AVAILABLE = False
    def njit(*args, **kwargs):
        if len(args) == 1 and callable(args[0]):
            return args[0]
        return lambda f: f


@njit(cache=True)
def _nn_residuals_jit(X, y, T, Z, dups, dupsid, matches, dT, dZ, n):
    """Inner loop of the NN-style variance helper.

    Pure-numeric kernel for `vce='nn'`. JIT-compiled when numba is
    installed; otherwise runs as-is in Python (slower but correct).
    All arrays are 1-D except Z which is (n, dZ); T and Z are dummy
    arrays of zeros when not used (caller guards with dT / dZ).
    """
    res = np.zeros((n, 1 + dT + dZ))
    for pos in range(n):
        rpos = dups[pos] - dupsid[pos]
        lpos = dupsid[pos] - 1
        cap = matches if matches < n - 1 else n - 1
        while lpos + rpos < cap:
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
        lo = pos - lpos
        hi = pos + rpos + 1
        Ji = (hi - lo) - 1
        sf = (Ji / (Ji + 1)) ** 0.5
        s_y = 0.0
        for j in range(lo, hi):
            s_y += y[j]
        s_y -= y[pos]
        res[pos, 0] = sf * (y[pos] - s_y / Ji)
        if dT == 1:
            s_T = 0.0
            for j in range(lo, hi):
                s_T += T[j]
            s_T -= T[pos]
            res[pos, 1] = sf * (T[pos] - s_T / Ji)
        for i in range(dZ):
            s_Z = 0.0
            for j in range(lo, hi):
                s_Z += Z[j, i]
            s_Z -= Z[pos, i]
            res[pos, 1 + dT + i] = sf * (Z[pos, i] - s_Z / Ji)
    return res


def resolve_name(arg, data, argname):
    """If `arg` is a column-name string and `data` is a DataFrame, return
    the corresponding column as a pandas Series. A list/tuple of names
    returns a DataFrame slice. Anything else passes through unchanged.
    """
    if arg is None:
        return None
    if isinstance(arg, str):
        if data is None:
            raise ValueError(
                f"String `{argname}` requires `data = ` (a pandas DataFrame)."
            )
        if arg not in data.columns:
            raise ValueError(
                f"Column {arg!r} not found in `data` (for `{argname}`)."
            )
        return data[arg]
    if (isinstance(arg, (list, tuple))
            and len(arg) > 0 and all(isinstance(c, str) for c in arg)):
        if data is None:
            raise ValueError(
                f"Character `{argname}` requires `data = ` (a pandas DataFrame)."
            )
        missing = [c for c in arg if c not in data.columns]
        if missing:
            raise ValueError(f"Column(s) not found in `data`: {missing}")
        return data[list(arg)]
    return arg


def resolve_covs(covs, data, argname='covs'):
    """Normalize the `covs` argument into a pandas DataFrame / 2-D array.

    Accepted forms:
      - None                                    -> None
      - str starting with '~' (formula)         -> model matrix via formulaic
                                                   (preferred) or patsy. The
                                                   intercept column, if
                                                   present, is dropped. NaN
                                                   rows are preserved so
                                                   downstream complete_cases
                                                   does the masking.
      - str or list/tuple of str (column names) -> resolve_name passthrough.
      - anything else (ndarray/DataFrame/etc.)  -> returned unchanged.
    """
    if covs is None:
        return None
    if isinstance(covs, (list, tuple)) and len(covs) == 0:
        raise ValueError(
            f"`{argname}` is empty. Supply at least one column name, "
            f"or omit `{argname}`."
        )
    if hasattr(covs, 'shape') and len(getattr(covs, 'shape', ())) >= 1 \
            and 0 in covs.shape:
        raise ValueError(
            f"`{argname}` has zero columns/rows. Supply at least one "
            f"covariate, or omit `{argname}`."
        )
    if isinstance(covs, str) and covs.strip().startswith('~'):
        if data is None:
            raise ValueError(
                f"Formula `{argname}` requires `data = ` (a pandas DataFrame)."
            )
        try:
            from formulaic import Formula
            mm = Formula(covs).get_model_matrix(data, na_action='ignore')
            df = mm if isinstance(mm, pd.DataFrame) else pd.DataFrame(mm)
        except ImportError:
            try:
                from patsy import dmatrix, NAAction
                df = dmatrix(covs, data, return_type='dataframe',
                             NA_action=NAAction(NA_types=[]))
            except ImportError:
                raise ImportError(
                    f"Formula `{argname}` requires `formulaic` or `patsy`. "
                    "Install via `pip install rdrobust[formula]` or "
                    "`pip install formulaic`."
                )
        drop = [c for c in df.columns if str(c).lower() in
                ('intercept', '(intercept)')]
        if drop:
            df = df.drop(columns=drop)
        return df
    return resolve_name(covs, data, argname)


def resolve_data_args(data, **kwargs):
    """Bulk-resolve a set of column-name args against `data`.

    Returns a dict mapping each kwarg name to the resolved value. If
    `data is None`, every value is passed through unchanged. Use as:

        loc = resolve_data_args(data, y=y, x=x, cluster=cluster, ...)
        y, x, cluster = loc['y'], loc['x'], loc['cluster']
    """
    if data is None:
        return kwargs
    return {name: resolve_name(val, data, name) for name, val in kwargs.items()}


def _capture_repr(print_func):
    buffer = io.StringIO()
    with redirect_stdout(buffer):
        print_func()
    return buffer.getvalue().rstrip()


class rdrobust_output:
    def __init__(self, Estimate, bws, coef, se, z, pv, ci, beta_p_l, beta_p_r,
                 V_cl_l, V_cl_r, V_rb_l, V_rb_r, N, N_h, N_b, M,
                 tau_cl, tau_bc, c, p, q, bias, kernel, all,
                 vce, bwselect, level, masspoints,
                 rdmodel=None, n_clust=None, n_clust_l=None, n_clust_r=None, coef_covs=None,
                 tau_T=None, se_T=None, z_T=None, pv_T=None, ci_T=None,
                 beta_T_p_l=None, beta_T_p_r=None):

        self.Estimate = Estimate
        self.bws = bws
        self.coef = coef
        self.se = se
        self.z = z
        self.pv = pv
        self.ci = ci
        self.beta_p_l = beta_p_l
        self.beta_p_r = beta_p_r
        self.V_cl_l = V_cl_l
        self.V_cl_r = V_cl_r
        self.V_rb_l = V_rb_l
        self.V_rb_r = V_rb_r
        self.N = N
        self.N_h = N_h
        self.N_b = N_b
        self.M = M
        self.tau_cl = tau_cl
        self.tau_bc = tau_bc
        self.c = c
        self.p = p
        self.q = q
        self.bias = bias
        self.kernel = kernel
        self.all = all
        self.vce = vce
        self.bwselect = bwselect
        self.level = level
        self.masspoints = masspoints
        self.rdmodel = rdmodel
        self.n_clust = n_clust
        self.n_clust_l = n_clust_l
        self.n_clust_r = n_clust_r
        self.coef_covs = coef_covs
        # First-stage (T) outputs for fuzzy designs; None for sharp designs.
        self.tau_T = tau_T
        self.se_T = se_T
        self.z_T = z_T
        self.pv_T = pv_T
        self.ci_T = ci_T
        self.beta_T_p_l = beta_T_p_l
        self.beta_T_p_r = beta_T_p_r

    def __repr__(self):
        return _capture_repr(self._repr_print)

    __str__ = __repr__

    def _repr_print(self):
        print('Call: rdrobust')
        fw = 30
        fw_r = 14

        # Model description
        if self.rdmodel is not None:
            print('')
            print(self.rdmodel)

        print('')
        print('Number of Observations:'.ljust(fw),
              str(self.N[0]+self.N[1]).rjust(fw_r))
        print('Polynomial Order Est. (p):'.ljust(fw), str(self.p).rjust(fw_r))
        print('Polynomial Order Bias (q):'.ljust(fw), str(self.q).rjust(fw_r))
        print('Kernel:'.ljust(fw), str(self.kernel).rjust(fw_r))
        print('Bandwidth Selection:'.ljust(fw), str(self.bwselect).rjust(fw_r))
        print('Var-Cov Estimator:'.ljust(fw), str(self.vce).rjust(fw_r))
        print('')

        # Table Sides
        fw_n = 25
        fw_l = 10
        fw_r = 10
        print(''.ljust(fw_n),'Left'.rjust(fw_l), 'Right'.rjust(fw_r))
        print('-'*(fw_n + fw_l + fw_r + 3))
        print('Number of Observations'.ljust(fw_n),
              str(self.N[0]).rjust(fw_l),
              str(self.N[1]).rjust(fw_r))
        print('Number of Unique Obs.'.ljust(fw_n),
              str(self.M[0]).rjust(fw_l),
              str(self.M[1]).rjust(fw_r))
        print('Number of Effective Obs.'.ljust(fw_n),
              str(self.N_h[0]).rjust(fw_l),
              str(self.N_h[1]).rjust(fw_r))
        print('Bandwidth Estimation'.ljust(fw_n),
              str(round(self.bws.iloc[0,0],3)).rjust(fw_l),
              str(round(self.bws.iloc[0,1],3)).rjust(fw_r))
        print('Bandwidth Bias'.ljust(fw_n),
              str(round(self.bws.iloc[1,0],3)).rjust(fw_l),
              str(round(self.bws.iloc[1,1],3)).rjust(fw_r))
        print('rho (h/b)'.ljust(fw_n),
              str(round(self.bws.iloc[0,0]/self.bws.iloc[1,0],3)).rjust(fw_l),
              str(round(self.bws.iloc[0,1]/self.bws.iloc[1,1],3)).rjust(fw_l))
        if self.n_clust is not None:
            # Per-side cluster counts (Left / Right)
            n_l = self.n_clust_l if self.n_clust_l is not None else self.n_clust
            n_r = self.n_clust_r if self.n_clust_r is not None else self.n_clust
            print('Number of Clusters'.ljust(fw_n),
                  str(n_l).rjust(fw_l),
                  str(n_r).rjust(fw_r))
        print('')

        # Table Estimation
        fw_l = 15
        fw_c = 8
        fw_ci = 18
        n_dec = 3
        print('Method'.ljust(fw_l),'Coef.'.rjust(fw_c),
              'S.E.'.rjust(fw_c),'z-stat'.rjust(fw_c),
              'P>|z|'.rjust(fw_c), (str(self.level) + '% CI').center(fw_ci))
        print('-'*73)

        if self.all is None: method_id = (0,2)
        else: method_id = (0,1,2)
        for j in method_id:
            if j==0:
                print(self.coef.index[j].ljust(fw_l),
                      str(round(self.coef.iloc[j].item(),n_dec)).rjust(fw_c),
                      str(round(self.se.iloc[j].item(),n_dec)).rjust(fw_c),
                      str(round(self.z.iloc[j].item(),n_dec)).rjust(fw_c),
                      ("{:.3e}".format(self.pv.iloc[j].item())).rjust(fw_c+3),
                      ('[' + str(round(self.ci.iloc[j,0].item(),n_dec)) + ', ' + str(round(self.ci.iloc[j,1].item(),n_dec)) + ']').rjust(fw_ci))
            else:
                if self.all:
                    print(self.coef.index[j].ljust(fw_l),
                      str(round(self.coef.iloc[j].item(),n_dec)).rjust(fw_c),
                      str(round(self.se.iloc[j].item(),n_dec)).rjust(fw_c),
                      str(round(self.z.iloc[j].item(),n_dec)).rjust(fw_c),
                      ("{:.3e}".format(self.pv.iloc[j].item())).rjust(fw_c+3),
                      ('[' + str(round(self.ci.iloc[j,0].item(),n_dec)) + ', ' + str(round(self.ci.iloc[j,1].item(),n_dec)) + ']').rjust(fw_ci))
                else:
                    print(self.coef.index[j].ljust(fw_l),
                      '-'.rjust(fw_c),
                      '-'.rjust(fw_c),
                      str(round(self.z.iloc[j].item(),n_dec)).rjust(fw_c),
                      ("{:.3e}".format(self.pv.iloc[j].item())).rjust(fw_c+3),
                      ('[' + str(round(self.ci.iloc[j,0].item(),n_dec)) + ', ' + str(round(self.ci.iloc[j,1].item(),n_dec)) + ']').rjust(fw_ci))

        # First-stage (T) panel for fuzzy designs, mirroring R's print.rdrobust.
        if self.tau_T is not None:
            print('')
            print('First-Stage Estimates.')
            print('Method'.ljust(fw_l),'Coef.'.rjust(fw_c),
                  'S.E.'.rjust(fw_c),'z-stat'.rjust(fw_c),
                  'P>|z|'.rjust(fw_c), (str(self.level) + '% CI').center(fw_ci))
            print('-'*73)
            for j in method_id:
                if j == 0:
                    print(self.tau_T.index[j].ljust(fw_l),
                      str(round(self.tau_T.iloc[j].item(), n_dec)).rjust(fw_c),
                      str(round(self.se_T.iloc[j].item(), n_dec)).rjust(fw_c),
                      str(round(self.z_T.iloc[j].item(),  n_dec)).rjust(fw_c),
                      ("{:.3e}".format(self.pv_T.iloc[j].item())).rjust(fw_c+3),
                      ('[' + str(round(self.ci_T.iloc[j,0].item(), n_dec)) + ', '
                           + str(round(self.ci_T.iloc[j,1].item(), n_dec)) + ']').rjust(fw_ci))
                else:
                    if self.all:
                        print(self.tau_T.index[j].ljust(fw_l),
                          str(round(self.tau_T.iloc[j].item(), n_dec)).rjust(fw_c),
                          str(round(self.se_T.iloc[j].item(), n_dec)).rjust(fw_c),
                          str(round(self.z_T.iloc[j].item(),  n_dec)).rjust(fw_c),
                          ("{:.3e}".format(self.pv_T.iloc[j].item())).rjust(fw_c+3),
                          ('[' + str(round(self.ci_T.iloc[j,0].item(), n_dec)) + ', '
                               + str(round(self.ci_T.iloc[j,1].item(), n_dec)) + ']').rjust(fw_ci))
                    else:
                        print(self.tau_T.index[j].ljust(fw_l),
                          '-'.rjust(fw_c),
                          '-'.rjust(fw_c),
                          str(round(self.z_T.iloc[j].item(),  n_dec)).rjust(fw_c),
                          ("{:.3e}".format(self.pv_T.iloc[j].item())).rjust(fw_c+3),
                          ('[' + str(round(self.ci_T.iloc[j,0].item(), n_dec)) + ', '
                               + str(round(self.ci_T.iloc[j,1].item(), n_dec)) + ']').rjust(fw_ci))

        return ''

class rdplot_output:
    def __init__(self, coef, rdplot, vars_bins, vars_poly, J, J_IMSE, J_MV,
                 scale, rscale, bin_avg, bin_med, p, c, h, N, N_h,
                 binselect, kernel, coef_covs=None):
        self.coef = coef
        self.rdplot = rdplot
        self.vars_bins = vars_bins
        self.vars_poly = vars_poly
        self.J = J
        self.J_IMSE = J_IMSE
        self.J_MV = J_MV
        self.scale = scale
        self.rscale = rscale
        self.bin_avg = bin_avg
        self.bin_med = bin_med
        self.p = p
        self.c = c
        self.h = h
        self.N = N
        self.N_h = N_h
        self.binselect = binselect
        self.kernel = kernel
        self.coef_covs = coef_covs
        
    def __repr__(self):
        return _capture_repr(self._repr_print)

    __str__ = __repr__

    def _repr_print(self):
        print('Call: rdplot')
        fw = 30
        fw_r = 14
        print('Number of Observations:'.ljust(fw),
              str(self.N[0]+self.N[1]).rjust(fw_r))
        print('Kernel:'.ljust(fw), str(self.kernel).rjust(fw_r))
        print('Polynomial Order Est. (p):'.ljust(fw), str(self.p).rjust(fw_r))
        print('')
        # Table Sides
        fw_n = 25
        fw_l = 10
        fw_r = 10
        nd = 3
        print(''.ljust(fw_n),'Left'.rjust(fw_l), 'Right'.rjust(fw_r))
        print('-'*(fw_n + fw_l + fw_r + 3))
        print('Number of Observations'.ljust(fw_n),
              str(self.N[0]).rjust(fw_l),
              str(self.N[1]).rjust(fw_r))
        print('Number of Effective Obs'.ljust(fw_n),
              str(self.N_h[0]).rjust(fw_l),
              str(self.N_h[1]).rjust(fw_r))
        print('Bandwidth poly. fit (h)'.ljust(fw_n),
              str(round(self.h[0],nd)).rjust(fw_l),
              str(round(self.h[1],nd)).rjust(fw_r))
        print('Number of bins scale'.ljust(fw_n),
              str(self.scale[0]).rjust(fw_l),
              str(self.scale[1]).rjust(fw_r))
        print('Bins Selected'.ljust(fw_n),
              str(self.J[0]).rjust(fw_l),
              str(self.J[1]).rjust(fw_r))
        print('Average Bin Length'.ljust(fw_n),
              str(round(self.bin_avg[0], nd)).rjust(fw_l),
              str(round(self.bin_avg[1], nd)).rjust(fw_r))
        print('Median Bin Length'.ljust(fw_n),
              str(round(self.bin_med[0], nd)).rjust(fw_l),
              str(round(self.bin_med[1], nd)).rjust(fw_r))
        print('IMSE-optimal bins'.ljust(fw_n),
              str(round(self.J_IMSE[0], nd)).rjust(fw_l),
              str(round(self.J_IMSE[1], nd)).rjust(fw_r))
        print('Mimicking Variance bins'.ljust(fw_n),
              str(round(self.J_MV[0], nd)).rjust(fw_l),
              str(round(self.J_MV[1], nd)).rjust(fw_r))
        print('')
        print('Relative to IMSE-optimal:'.ljust(fw_n))
        print('Implied scale'.ljust(fw_n),
              str(round(self.rscale[0], nd)).rjust(fw_l),
              str(round(self.rscale[1], nd)).rjust(fw_r))
        print('WIMSE variance weight'.ljust(fw_n),
              str(round(1/(1+self.rscale[0]**3),nd)).rjust(fw_l),
              str(round(1/(1+self.rscale[1]**3),nd)).rjust(fw_r))
        print('WIMSE bias weight'.ljust(fw_n),
              str(round(self.rscale[0]**3/(1+self.rscale[0]**3),nd)).rjust(fw_l),
              str(round(self.rscale[1]**3/(1+self.rscale[1]**3),nd)).rjust(fw_r))
        return ''

class rdbwselect_output:
    def __init__(self, bws, bwselect, kernel, p, q, c,
                 N, N_h, M, vce, masspoints):
        self.bws = bws
        self.bwselect = bwselect
        self.kernel = kernel
        self.p = p
        self.q = q
        self.c = c
        self.N = N
        self.N_h = N_h
        self.M = M
        self.vce = vce
        self.masspoints = masspoints
    
    def __repr__(self):
        return _capture_repr(self._repr_print)

    __str__ = __repr__

    def _repr_print(self, nd = 3):
        fw = 30
        fw_r = 14
        print('Call: rdbwselect')
        print('Number of Observations:'.ljust(fw),
              str(self.N[0]+self.N[1]).rjust(fw_r))
        print('Polynomial Order Est. (p):'.ljust(fw), str(self.p).rjust(fw_r))
        print('Polynomial Order Bias (q):'.ljust(fw), str(self.q).rjust(fw_r))
        print('Kernel:'.ljust(fw), str(self.kernel).rjust(fw_r))
        print('Bandwidth Selection:'.ljust(fw), str(self.bwselect).rjust(fw_r))
        print('Var-Cov Estimator:'.ljust(fw), str(self.vce).rjust(fw_r))
        print('')
        print(self.bws.round(nd))
        return ''
        
def tomat(x):
    return x.reshape(len(x),-1);

def ncol(x):
    try: d = x.shape[1]
    except: d = 1
    return d;

def crossprod(x, y = None):
    if y is None: return np.matmul(x.T,x)
    else: return np.matmul(x.T,y)

def quantile_type2(x, probs):
    """Hyndman-Fan type 2 quantiles, matching R's quantile(..., type=2)."""
    x = np.sort(np.asarray(x, dtype=float).reshape(-1))
    n = x.size
    if n == 0:
        raise ValueError("cannot compute quantiles of an empty array")

    probs_arr = np.asarray(probs, dtype=float)
    scalar = probs_arr.ndim == 0
    probs_flat = probs_arr.reshape(-1)
    out = np.empty(probs_flat.size, dtype=float)

    for idx, prob in enumerate(probs_flat):
        if prob <= 0:
            out[idx] = x[0]
        elif prob >= 1:
            out[idx] = x[-1]
        else:
            h = n * prob
            j = int(np.floor(h))
            if h == j:
                out[idx] = 0.5 * (x[j - 1] + x[j])
            else:
                out[idx] = x[j]

    if scalar:
        return out[0]
    return out.reshape(probs_arr.shape)

def _vander(u, p):
    """Q1: build Vandermonde [1, u, u^2, ..., u^p] via successive
    multiplication. Avoids np.power on every column and reuses memory."""
    u = np.asarray(u).reshape(-1)
    n = u.size
    out = np.empty((n, p + 1))
    out[:, 0] = 1.0
    if p >= 1:
        for j in range(1, p + 1):
            out[:, j] = out[:, j - 1] * u
    return out


def _build_cluster_idx(C_flat):
    """Build the per-cluster row-index list once. T3 helper: callers
    that invoke rdrobust_vce multiple times with the same C should cache
    the result of this and pass it via cluster_idx=."""
    clusters = np.unique(C_flat)
    sort_idx = np.argsort(C_flat)
    C_sorted = C_flat[sort_idx]
    split_pts = np.searchsorted(C_sorted, clusters, side='right')
    start_pts = np.concatenate([[0], split_pts[:-1]])
    return [sort_idx[start_pts[i]:split_pts[i]] for i in range(len(clusters))]

def nanmat(n, m = None):
    # Create a (n x m) matrix of NaN
    if m is None: M = np.empty((n,))
    else: M = np.empty((n,m,))    
    M.fill(np.nan)
    return M;

def inv_chol(x):
    # Inverse of a symmetric positive-definite matrix via Cholesky + triangular solve.
    # No check is made if x is indeed positive definite!
    c, low = cho_factor(x, lower=True)
    return cho_solve((c, low), np.eye(x.shape[0]))
        
def qrXXinv(x):
    return inv_chol(crossprod(x,x))

def complete_cases(x):
    return np.all(np.invert(np.isnan(x)), axis = 1)

def covs_drop_fun(z, tol=1e-5):
    """Drop collinear columns from `z`.

    Fast path (the common case): if `z` is full column rank, return it
    unchanged. We detect this via numpy's QR (orders of magnitude faster
    than scipy's column-pivoted QR for tall-thin matrices), checking for
    near-zero diagonal entries in R.

    Slow path (rank-deficient z): fall back to scipy's column-pivoted QR
    so we can identify which specific columns to drop.
    """
    q, r = np.linalg.qr(z)
    diag_abs = np.abs(np.diagonal(r))
    if np.all(diag_abs > tol):
        return z
    q, r, pivot = qr(a=z, pivoting=True)
    keep = pivot[np.abs(np.diagonal(r)) > tol]
    return z[:, keep]
    
def rdrobust_kweight(X, c, h, kernel):
    
    u = (X-c)/h
    if kernel=="epanechnikov" or kernel=="epa":
        w = (0.75*(1-u**2)*(abs(u)<=1))/h
    elif kernel=="uniform" or kernel=="uni":
         w = (0.5*(abs(u)<=1))/h
    else:
        w = ((1-abs(u))*(abs(u)<=1))/h
    return w

def rdrobust_res(X, y, T, Z, m, hii, vce, matches, dups, dupsid, d,
                 crv3=False, crv2=False, has_cluster=False):

    n = len(y)
    dT = dZ =  0
    if T is not None:
        dT = 1
    if Z is not None:
        dZ = ncol(Z)
    res = nanmat(n,1+dT+dZ)
    if vce=="nn":
        X_arr = np.ascontiguousarray(X, dtype=np.float64).ravel()
        y_arr = np.ascontiguousarray(y, dtype=np.float64).ravel()
        T_arr = (np.ascontiguousarray(T, dtype=np.float64).ravel()
                 if T is not None else np.zeros(n, dtype=np.float64))
        Z_arr = (np.ascontiguousarray(Z, dtype=np.float64).reshape(n, dZ)
                 if Z is not None else np.zeros((n, 0), dtype=np.float64))
        dups_arr   = np.ascontiguousarray(dups,   dtype=np.int64).ravel()
        dupsid_arr = np.ascontiguousarray(dupsid, dtype=np.int64).ravel()
        res = _nn_residuals_jit(X_arr, y_arr, T_arr, Z_arr,
                                dups_arr, dupsid_arr,
                                int(matches), int(dT), int(dZ), int(n))
    elif crv3 or crv2:
        # CRV2/CRV3 mode: return raw (unweighted) residuals.
        # Cluster-level hat-matrix adjustment is handled in rdrobust_vce.
        m = m.reshape(-1,ncol(m))
        res[:,0] = y - m[:,0]
        if dT==1: res[:,1] = T - m[:,1]
        if dZ>0:
            for i in range(dZ):
                res[:,1+dT+i] = Z[:,i] - m[:,1+dT+i]
    else:

        if vce=="hc0": w = 1
        elif vce=="hc1": w = 1 if has_cluster else np.sqrt(n/(n-d))
        elif vce=="hc2":
            hii = hii.reshape(-1)
            w = np.sqrt(1/np.maximum(1-hii, 1e-8))
        else:
            hii = hii.reshape(-1)
            w = 1/np.maximum(1-hii, 1e-8)
        m = m.reshape(-1,ncol(m))
        res[:,0] = w*(y-m[:,0])
        if dT==1: res[:,1] = w*(T-m[:,1])
        if dZ>0:
           for i in range(dZ):
               res[:,1+dT+i] = w*(Z[:,i]- m[:,1+dT+i])

    return res

def rdrobust_vce(d, s, RX, res, C, invG=None, sqrtRX=None, crv2=False,
                 k_override=None, cluster_idx=None, G=None):
    # G: optional precomputed X'WX (the Gram). When supplied, CRV3 skips the
    # G = inv(invG) round-trip the caller already had to compute X'WX before
    # forming invG, so passing it through avoids two unnecessary k×k inverses
    # per call plus the inv(inv()) numerical round-trip error.
    # k_override: when supplied, overrides the (n-1)/(n-k) df correction's k
    # for CR1. Used when RX is a "score-like" matrix (e.g. Q_q) whose ncol is
    # smaller than the effective number of regressors estimated. Without it,
    # CR1 at h!=b under cluster used k=p+1 from ncol(Q_q), inconsistent with
    # the q-regression direct path at h=b which uses k=q+1.
    # cluster_idx: optional precomputed (sort_idx, start_pts, split_pts).
    # When supplied, skips the per-call np.argsort/np.searchsorted on C.
    # Used by callers that invoke rdrobust_vce multiple times with the same C.
    # Note: k is the matrix dimension; k_df is what enters the df correction.
    k    = ncol(RX)
    k_df = k if k_override is None else k_override
    M = np.zeros((k,k))
    crv3_mode = invG is not None and C is not None and not crv2
    crv2_mode = invG is not None and C is not None and crv2

    if C is None:
        w = 1
        if d==0:
            M  = crossprod(res*RX,res*RX)
        else:
            # T4: Σᵢⱼ sᵢsⱼ rᵢrⱼ factors as (Σₗ sₗ rₗ)². One weighted
            # crossprod instead of an (d+1)² Python loop. Algebra: with
            # r_comb[n] = Σₗ s[l]*res[n,l], we have
            #   sum_{i,j} s[i] s[j] res[:,i] res[:,j] = r_comb**2,
            # and crossprod(RX * r_comb.reshape(-1,1)) produces
            # sum_n r_comb[n]² · RX[n,:]ᵀ RX[n,:], which is the M we want.
            s_arr  = np.asarray(s).reshape(-1)
            r_comb = res @ s_arr
            M = crossprod(RX * r_comb.reshape(-1, 1))
    elif crv3_mode:
        # CRV3: cluster-level hat-matrix adjustment (Pustejovsky & Tipton 2018).
        C_flat = C.reshape(-1)
        # T3: use cached cluster_idx if supplied, else build once.
        cl_indices = _build_cluster_idx(C_flat) if cluster_idx is None else cluster_idx
        g = len(cl_indices)
        w = 1  # CRV3 is approximately unbiased; no df correction needed
        has_sqrtRX = sqrtRX is not None
        # Py-2: prefer caller-supplied G (= X'WX) to avoid inv(inv(.)) round-trip.
        if G is None:
            try:
                G = np.linalg.inv(invG)
            except:
                G = np.linalg.pinv(invG)

        # Py-3/Py-4: M_g = inv(G - L_g) only ever consumed as M_g @ u_g_l, so use
        # solve (d=0: single use) or lu_factor + lu_solve (d>0: reuse across l-loop)
        # instead of forming the explicit inverse.
        if d==0:
            for i in range(g):
                ind = cl_indices[i]
                if has_sqrtRX:
                    sRX_g = sqrtRX[ind,:]
                    sw_g  = sRX_g[:,0:1]
                    e_g   = res[ind,0:1]
                    L_g   = crossprod(sRX_g)
                    u_g   = crossprod(sRX_g, sw_g * e_g)
                    try:
                        adj = np.linalg.solve(G - L_g, u_g)
                    except:
                        adj = np.linalg.pinv(G - L_g) @ u_g
                    score_g = u_g + L_g @ adj
                else:
                    RX_g = RX[ind,:]
                    e_g  = res[ind,0:1]
                    L_g  = crossprod(RX_g)
                    u_g  = crossprod(RX_g, e_g)
                    try:
                        adj = np.linalg.solve(G - L_g, u_g)
                    except:
                        adj = np.zeros_like(u_g)
                    score_g = u_g + L_g @ adj
                M = M + score_g @ score_g.T
        else:
            # Py-4 note: for very small k (typical RD: 4-8), scipy's
            # lu_factor/lu_solve Python overhead exceeds the flop savings vs
            # np.linalg.inv on covs_3+CR3. Keep the explicit inverse here;
            # it's reused d+1 times across the residual columns anyway.
            for i in range(g):
                ind = cl_indices[i]
                ri  = res[ind,:]
                sv  = np.zeros(k)
                if has_sqrtRX:
                    sRX_g = sqrtRX[ind,:]
                    sw_g  = sRX_g[:,0:1]
                    L_g   = crossprod(sRX_g)
                    try:
                        M_g = np.linalg.inv(G - L_g)
                    except:
                        M_g = np.linalg.pinv(G - L_g)
                    for l in range(d+1):
                        u_g_l = crossprod(sRX_g, sw_g * ri[:,l:l+1])
                        score_l = u_g_l + L_g @ (M_g @ u_g_l)
                        sv      = sv + score_l.reshape(-1) * s[l]
                else:
                    RX_g = RX[ind,:]
                    L_g  = crossprod(RX_g)
                    try:
                        M_g = np.linalg.inv(G - L_g)
                    except:
                        M_g = np.zeros((k,k))
                    for l in range(d+1):
                        u_g_l = crossprod(RX_g, ri[:,l:l+1])
                        score_l = u_g_l + L_g @ (M_g @ u_g_l)
                        sv      = sv + score_l.reshape(-1) * s[l]
                M = M + np.outer(sv, sv)
    elif crv2_mode:
        # CRV2: Bell & McCaffrey (2002) leverage-adjusted cluster-robust.
        C_flat = C.reshape(-1)
        # T3: use cached cluster_idx if supplied, else build once.
        cl_indices = _build_cluster_idx(C_flat) if cluster_idx is None else cluster_idx
        g = len(cl_indices)
        w = 1  # CRV2 is approximately unbiased
        has_sqrtRX = sqrtRX is not None
        try:
            C_half  = np.linalg.cholesky(invG).T  # upper Cholesky
            tC_half = C_half.T
            crv2_ok = True
        except:
            crv2_ok = False

        if d==0:
            for i in range(g):
                ind = cl_indices[i]
                if has_sqrtRX:
                    sRX_g = sqrtRX[ind,:]
                    sw_g  = sRX_g[:,0:1]
                    e_g   = res[ind,0:1]
                    L_g   = crossprod(sRX_g)
                    u_g   = crossprod(sRX_g, sw_g * e_g)
                else:
                    RX_g = RX[ind,:]
                    e_g  = res[ind,0:1]
                    L_g  = crossprod(RX_g)
                    u_g  = crossprod(RX_g, e_g)
                if crv2_ok:
                    F_sq   = C_half @ L_g @ tC_half
                    sigma2, V_g = np.linalg.eigh(F_sq)
                    sigma2 = np.maximum(sigma2, 0)
                    c_coef = np.where(sigma2 < 1e-14, 0,
                                      (1/np.sqrt(np.maximum(1-sigma2, 1e-8)) - 1) / sigma2)
                    Ru_g    = C_half @ u_g
                    adj     = V_g @ (c_coef.reshape(-1,1) * (V_g.T @ Ru_g))
                    score_g = u_g + L_g @ (tC_half @ adj)
                else:
                    score_g = u_g
                M = M + score_g @ score_g.T
        else:
            for i in range(g):
                ind = cl_indices[i]
                ri  = res[ind,:]
                sv  = np.zeros(k)
                if has_sqrtRX:
                    sRX_g = sqrtRX[ind,:]
                    sw_g  = sRX_g[:,0:1]
                    L_g   = crossprod(sRX_g)
                else:
                    RX_g = RX[ind,:]
                    L_g  = crossprod(RX_g)
                if crv2_ok:
                    F_sq     = C_half @ L_g @ tC_half
                    sigma2, V_g = np.linalg.eigh(F_sq)
                    sigma2   = np.maximum(sigma2, 0)
                    c_coef   = np.where(sigma2 < 1e-14, 0,
                                        (1/np.sqrt(np.maximum(1-sigma2, 1e-8)) - 1) / sigma2)
                    LtC_half = L_g @ tC_half
                    for l in range(d+1):
                        if has_sqrtRX:
                            u_g_l = crossprod(sRX_g, sw_g * ri[:,l:l+1])
                        else:
                            u_g_l = crossprod(RX_g, ri[:,l:l+1])
                        Ru_g_l  = C_half @ u_g_l
                        adj_l   = V_g @ (c_coef.reshape(-1,1) * (V_g.T @ Ru_g_l))
                        score_l = u_g_l + LtC_half @ adj_l
                        sv      = sv + score_l.reshape(-1) * s[l]
                else:
                    for l in range(d+1):
                        if has_sqrtRX:
                            u_g_l = crossprod(sRX_g, sw_g * ri[:,l:l+1])
                        else:
                            u_g_l = crossprod(RX_g, ri[:,l:l+1])
                        sv = sv + u_g_l.reshape(-1) * s[l]
                M = M + np.outer(sv, sv)
    else:
        # Standard cluster-robust (CR1) with small-sample df correction
        C_flat = C.reshape(-1)
        try: n = len(C_flat)
        except: n = 1
        # T3: use cached cluster_idx if supplied, else build once.
        cl_indices = _build_cluster_idx(C_flat) if cluster_idx is None else cluster_idx
        g = len(cl_indices)
        w = ((n-1)/(n-k_df))*(g/(g-1))

        if d==0:
            for i in range(g):
                ind = cl_indices[i]
                Xi = RX[ind,:]
                ri = res[ind,:]
                Xr = crossprod(Xi,ri).T
                M = M + crossprod(Xr,Xr)
        else:
            for i in range(g):
                ind = cl_indices[i]
                Xi = RX[ind,:]
                ri = res[ind,:]
                MHolder = np.zeros((1+d,k))
                for l in range(d+1):
                    MHolder[l,:] = crossprod(Xi,s[l]*ri[:,l]).T
                summedvalues = np.sum(MHolder, axis = 0)
                M = M + np.outer(summedvalues,summedvalues)

    return w*M


def rdrobust_vce_qq_cluster(Q, R_q, W_b, invG_q, res, C,
                            crv2=False, s=1, d=0, cluster_idx=None,
                            G_q=None):
    """Decoupled cluster-robust variance for the bias-corrected estimator when h != b.

    Sandwich uses Q_q (mixed-bandwidth) but leverage comes from the q-regression
    (bandwidth b), whose cluster hat matrix has eigvals in [0, 1). Matrix analog
    of using hii_q for HC2/HC3 in the same h != b branch.

    cluster_idx: optional precomputed _build_cluster_idx(C_flat) — when
    multiple calls share C, pass it through to avoid the O(N log N) sort.

    Returns meat matrix M such that V_rb = invG_p @ M @ invG_p.
    """
    Q      = np.asarray(Q, dtype=float)
    R_q    = np.asarray(R_q, dtype=float)
    W_b    = np.asarray(W_b, dtype=float).reshape(-1)
    invG_q = np.asarray(invG_q, dtype=float)
    res    = np.asarray(res, dtype=float)
    if res.ndim == 1:
        res = res.reshape(-1, 1)
    k   = Q.shape[1]
    k_R = R_q.shape[1]
    M   = np.zeros((k, k))

    C_flat = np.asarray(C).reshape(-1)
    # T3: use cached cluster_idx if supplied, else build once.
    cl_indices = _build_cluster_idx(C_flat) if cluster_idx is None else cluster_idx
    g = len(cl_indices)

    # Py-2: prefer caller-supplied G_q (= R_q'W_b R_q) to skip inv(inv(.)) round-trip.
    if G_q is None:
        try:
            G_q = np.linalg.inv(invG_q)
        except np.linalg.LinAlgError:
            G_q = np.linalg.pinv(invG_q)
    else:
        G_q = np.asarray(G_q, dtype=float)

    # s normalized to array for vectorized access
    if np.isscalar(s):
        s_arr = np.array([s], dtype=float)
    else:
        s_arr = np.asarray(s, dtype=float).reshape(-1)

    for gi in range(g):
        idx  = cl_indices[gi]
        Q_g  = Q[idx, :]
        R_g  = R_q[idx, :]
        Wb_g = W_b[idx]
        sR_g = R_g * np.sqrt(Wb_g).reshape(-1, 1)        # n_g x k_R, transformed
        L_g  = sR_g.T @ sR_g                              # k_R x k_R
        P_g  = Q_g.T @ R_g                                # k x k_R

        # T2: per-cluster Lambda (CRV2) or M_g (CRV3) depends only on
        # L_g / G_q / invG_q -- not on the residual column l. Compute once
        # per cluster, reuse across the l-loop. (Py-4: lu_factor/lu_solve was
        # tried but lost to np.linalg.inv on the typical k ≤ 8 RD design due
        # to per-call scipy overhead, so we keep the explicit inverse here.)
        Lambda = None
        M_g    = None
        use_simple = False
        if crv2:
            try:
                sigma2, V_sv = np.linalg.eigh(L_g)
            except np.linalg.LinAlgError:
                use_simple = True
            else:
                sigma2 = np.maximum(sigma2, 0.0)
                sigma  = np.sqrt(sigma2)
                tol = (sigma.max() * 1e-10) if sigma.size else 0.0
                keep = sigma > tol
                if not np.any(keep):
                    use_simple = True
                else:
                    V_r = V_sv[:, keep]
                    sig = sigma[keep]
                    Dsig = np.diag(sig)
                    T_mat = Dsig @ (V_r.T @ (invG_q @ V_r)) @ Dsig
                    T_mat = 0.5 * (T_mat + T_mat.T)
                    tau, W_g = np.linalg.eigh(T_mat)
                    tau = np.minimum(np.maximum(tau, 0.0), 1 - 1e-10)
                    gamma = 1.0 / np.sqrt(1.0 - tau) - 1.0
                    Dinv = np.diag(1.0 / sig)
                    A = V_r @ Dinv @ W_g
                    Lambda = A @ (gamma.reshape(-1, 1) * A.T)
        else:
            try:
                M_g = np.linalg.inv(G_q - L_g)
            except np.linalg.LinAlgError:
                M_g = np.zeros((k_R, k_R))

        sv = np.zeros(k)
        for l in range(d + 1):
            e_gl = res[idx, l]
            u_Qq = Q_g.T @ e_gl
            if use_simple:
                score_l = u_Qq.copy()
            else:
                u_q  = R_g.T @ (Wb_g * e_gl)
                adj  = Lambda @ u_q if crv2 else M_g @ u_q
                score_l = u_Qq + P_g @ adj
            sv = sv + score_l * s_arr[l]

        M = M + np.outer(sv, sv)

    return M

def _rdrobust_bw_Vfit(Y, X, T, Z, C, W, c, o, nu, h_V,
                      vce, nnmatch, kernel, dups, dupsid, covs_drop_coll):
    """V-fit portion of rdrobust_bw. Returns (V_V, BConst, s).

    Depends only on (side-data, o, nu, h_V, vce). In rdbwselect h_V is a fixed
    pilot bandwidth shared across calls, so results can be cached by (o, nu).
    """
    crv3 = (vce=="crv3") and C is not None
    crv2 = (vce=="crv2") and C is not None

    dT = dZ = 0
    w = rdrobust_kweight(X, c, h_V, kernel)
    if not np.isscalar(W): w = W*w
    ind_V = w> 0
    eY = Y[ind_V]
    eX = X[ind_V]
    eW = w[ind_V]
    n_V = np.sum(ind_V)
    D_V = eY.copy()
    R_V = _vander(eX - c, o)
    # Py-2: compute G = R'WR once, then invG via inv_chol. Lets us pass G to
    # rdrobust_vce so the CRV3 branch can skip the inv(invG) round-trip.
    RWsq_V = R_V * np.sqrt(eW).reshape(-1,1)
    G_V    = crossprod(RWsq_V)
    invG_V = inv_chol(G_V)
    s = 1
    eT = eC = eZ = None
    if T is not None:
        dT = 1
        eT = T[ind_V]
        D_V = np.column_stack((D_V,eT))
    if Z is not None:
        dZ = ncol(Z)
        eZ = Z[ind_V,:]
        D_V = np.column_stack((D_V,eZ))
        U = crossprod(R_V*eW.reshape(-1,1),D_V)
        ZWD  = crossprod(eZ*eW.reshape(-1,1),D_V)
        colsZ = np.arange(1+dT, max(2+dT+dZ-1,2+dT))
        UiGU =  crossprod(U[:,colsZ],np.matmul(invG_V,U))
        ZWZ = ZWD[:,colsZ] - UiGU[:,colsZ]
        ZWY = ZWD[:,:(1+dT)] - UiGU[:,:(1+dT)]
        if covs_drop_coll == 1:
            gamma = np.matmul(np.linalg.pinv(ZWZ),ZWY)
        else:
            gamma = np.matmul(inv_chol(ZWZ),ZWY)
        s = np.append(1, -gamma[:,0])

    if C is not None: eC =  C[ind_V]

    beta_V = np.matmul(invG_V,crossprod(R_V*eW.reshape(-1,1),D_V))
    if Z is None and T is not None:
        tau_Y = math.factorial(nu)*beta_V[nu,0]
        tau_T = math.factorial(nu)*beta_V[nu,1]
        s = np.array([1/tau_T , -(tau_Y/tau_T**2)])

    if Z is not None and T is not None:
        s_T = np.append(1, -gamma[:,1])
        beta_Y = np.append(beta_V[nu,0],beta_V[nu,colsZ])
        tau_Y = math.factorial(nu)*np.matmul(s, beta_Y)
        beta_T = np.append(beta_V[nu,1],beta_V[nu,colsZ])
        tau_T = math.factorial(nu)*np.matmul(s_T, beta_T)
        s = (np.hstack([1/tau_T , -(tau_Y/tau_T**2) ,
            -(1/tau_T)*gamma[:,0] + (tau_Y/tau_T**2)*gamma[:,1]]))

    dups_V = dupsid_V = hii = predicts_V = 0
    if vce == "nn":
        dups_V   = dups[ind_V]
        dupsid_V = dupsid[ind_V]

    if (vce in ("hc0","hc1","hc2","hc3")) and not crv3 and not crv2:
        predicts_V = np.matmul(R_V,beta_V)
        if vce=="hc2" or vce=="hc3":
            hii = np.sum(np.matmul(R_V,invG_V)*(R_V*eW.reshape(-1,1)), axis = 1).reshape(-1,1)
    elif crv3 or crv2:
        predicts_V = np.matmul(R_V,beta_V)

    sqrtRX_V  = RWsq_V if (crv3 or crv2) else None
    invG_V_c  = invG_V if (crv3 or crv2) else None
    G_V_c     = G_V    if crv3 else None  # Py-2: CRV2 only needs invG (via cholesky)
    res_V = (rdrobust_res(eX, eY, eT, eZ, predicts_V, hii, vce,
                          nnmatch, dups_V, dupsid_V, o+1,
                          crv3=crv3, crv2=crv2, has_cluster=C is not None))

    aux = rdrobust_vce(dT+dZ, s, R_V*eW.reshape(-1,1), res_V, eC,
                       invG=invG_V_c, sqrtRX=sqrtRX_V, crv2=crv2, G=G_V_c)
    V_V = np.matmul(np.matmul(invG_V,aux),invG_V)[nu,nu]

    v = crossprod(R_V*eW.reshape(-1,1),((eX-c)/h_V)**(o+1))
    Hp = np.zeros(o+1)
    for j in range(o+1): Hp[j] = h_V**j
    BConst = (Hp*np.matmul(invG_V,v))[nu]

    return V_V, BConst, s


def rdrobust_bw(Y, X, T, Z, C, W, c, o, nu, o_B, h_V, h_B, scale,
                 vce, nnmatch, kernel, dups, dupsid, covs_drop_coll,
                 _vcache=None):

    crv3 = (vce=="crv3") and C is not None
    crv2 = (vce=="crv2") and C is not None

    dT = 1 if T is not None else 0
    dZ = ncol(Z) if Z is not None else 0

    key = (o, nu)
    if _vcache is not None and key in _vcache:
        V_V, BConst, s = _vcache[key]
    else:
        V_V, BConst, s = _rdrobust_bw_Vfit(Y, X, T, Z, C, W, c, o, nu, h_V,
                                           vce, nnmatch, kernel, dups, dupsid,
                                           covs_drop_coll)
        if _vcache is not None:
            _vcache[key] = (V_V, BConst, s)

    w = rdrobust_kweight(X, c, h_B, kernel)
    if not np.isscalar(W): w = W*w
    ind = w> 0
    n_B = sum(ind)
    eY = Y[ind]
    eX = X[ind]
    eW = w[ind]
    D_B = eY
    R_B = _vander(eX - c, o_B)
    # Py-2: G_B kept around for the CRV3 branch.
    RWsq_B = R_B * np.sqrt(eW).reshape(-1,1)
    G_B    = crossprod(RWsq_B)
    invG_B = inv_chol(G_B)

    eT = eC = eZ = None
    if T is not None:
        eT = T[ind]
        D_B = np.column_stack((D_B,eT))
    if Z is not None:
        eZ = Z[ind,:]
        D_B = np.column_stack((D_B,eZ))
    if C is not None:
        eC = C[ind]

    beta_B = np.matmul(invG_B,crossprod(R_B*eW.reshape(-1,1),D_B))
    BWreg=0
    if scale>0:
        e_B = np.zeros(o_B+1)
        e_B[o+1]=1
        dups_B = dupsid_B = hii = predicts_B = 0
        if vce=="nn":
            dups_B   = dups[ind]
            dupsid_B = dupsid[ind]
        if (vce in ("hc0","hc1","hc2","hc3")) and not crv3 and not crv2:
            predicts_B = np.matmul(R_B,beta_B)
            if vce=="hc2" or vce=="hc3":
                hii = np.sum(np.matmul(R_B,invG_B)*(R_B*eW.reshape(-1,1)), axis = 1).reshape(-1,1)
        elif crv3 or crv2:
            predicts_B = np.matmul(R_B,beta_B)

        sqrtRX_B = RWsq_B if (crv3 or crv2) else None
        invG_B_c = invG_B if (crv3 or crv2) else None
        G_B_c    = G_B    if crv3 else None  # Py-2
        res_B = (rdrobust_res(eX, eY, eT, eZ, predicts_B, hii, vce, nnmatch,
                              dups_B, dupsid_B, o_B+1,
                              crv3=crv3, crv2=crv2, has_cluster=C is not None))
        aux = rdrobust_vce(dT+dZ, s, R_B*eW.reshape(-1,1), res_B, eC,
                           invG=invG_B_c, sqrtRX=sqrtRX_B, crv2=crv2, G=G_B_c)
        V_B = np.matmul(np.matmul(invG_B,aux),invG_B)[-1,-1]
        BWreg = 3*BConst**2*V_B
    try: beta_aux = beta_B[-1,:]
    except: beta_aux = beta_B[-1]
    B =  np.sqrt(2*(o+1-nu))*BConst*np.dot(s,beta_aux)
    V = (2*nu+1)*h_V**(2*nu+1)*V_V
    R = scale*(2*(o+1-nu))*BWreg
    rate = 1/(2*o+3)

    return V, B, R, rate
