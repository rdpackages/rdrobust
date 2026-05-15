#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 19:01:48 2021

@author: rmasini
"""

import numpy as np
import pandas  as pd
import math
import scipy.stats as sct
from rdrobust.rdbwselect import rdbwselect, _rdbwselect_compute
from rdrobust.funs import *
     
def rdrobust(y, x, c = None, fuzzy = None, deriv = None,
             p = None, q = None, h = None, b = None, rho = None,
             covs = None, covs_drop = True,
             kernel = "tri", weights = None, bwselect = "mserd",
             vce = "nn", cluster = None, nnmatch = 3, level = 95,
             scalepar = 1, scaleregul = 1, sharpbw = False,
             all = None, subset = None, masspoints = "adjust",
             bwcheck = None, bwrestrict = True, stdvars = False,
             data = None):
    
    '''
    Implements local polynomial Regression Discontinuity (RD) point estimators with robust bias-corrected confidence intervals and inference procedures developed in Calonico, Cattaneo and Titiunik (2014a), Calonico, Cattaneo and Farrell (2018), Calonico, Cattaneo, Farrell and Titiunik (2019), and Calonico, Cattaneo and Farrell (2020). It also computes alternative estimation and inference procedures available in the literature.
   
    Parameters
    ----------
    y: array	
    is the dependent variable.
    
    x	
    is the running variable (a.k.a. score or forcing variable).
    
    c	
    specifies the RD cutoff in x; default is c = 0.
    
    fuzzy	
    specifies the treatment status variable used to implement fuzzy RD estimation (or Fuzzy Kink RD if deriv=1 is also specified). Default is Sharp RD design and hence this option is not used.
    
    deriv	
    specifies the order of the derivative of the regression functions to be estimated. Default is deriv=0 (for Sharp RD, or for Fuzzy RD if fuzzy is also specified). Setting deriv=1 results in estimation of a Kink RD design (up to scale), or Fuzzy Kink RD if fuzzy is also specified.
    
    p	
    specifies the order of the local-polynomial used to construct the point-estimator; default is p = 1 (local linear regression).
    
    q	
    specifies the order of the local-polynomial used to construct the bias-correction; default is q = 2 (local quadratic regression).
    
    h	
    specifies the main bandwidth used to construct the RD point estimator. If not specified, bandwidth h is computed by the companion command rdbwselect. If two bandwidths are specified, the first bandwidth is used for the data below the cutoff and the second bandwidth is used for the data above the cutoff.
    
    b	
    specifies the bias bandwidth used to construct the bias-correction estimator. If not specified, bandwidth b is computed by the companion command rdbwselect. If two bandwidths are specified, the first bandwidth is used for the data below the cutoff and the second bandwidth is used for the data above the cutoff.
    
    rho	
    specifies the value of rho, so that the bias bandwidth b equals h/rho. Default is rho = 1 if h is specified but b is not.
    
    covs
    specifies additional covariates to be used for estimation and inference.
    Accepted forms:
      - a numpy array, pandas DataFrame, or Series (back-compatible);
      - a one-sided formula string (requires `data`), e.g.
        "~ z1 + z2 + I(z3**2)": expanded via `formulaic` (preferred) or
        `patsy`; the intercept column is dropped automatically;
      - a column name or list of column names (requires `data`),
        e.g. ["z1", "z2"]: resolved as `data[covs]`.

    covs_drop
    if TRUE, it checks for collinear additional covariates and drops them. Default is TRUE.
    
    kernel	
    is the kernel function used to construct the local-polynomial estimator(s). Options are triangular (default option), epanechnikov and uniform.
    
    weights	
    is the variable used for optional weighting of the estimation procedure. The unit-specific weights multiply the kernel function.
    
    bwselect	
    specifies the bandwidth selection procedure to be used. By default it computes both h and b, unless rho is specified, in which case it only computes h and sets b=h/rho.
    
    Options are:
    
    mserd one common MSE-optimal bandwidth selector for the RD treatment effect estimator.
    
    msetwo two different MSE-optimal bandwidth selectors (below and above the cutoff) for the RD treatment effect estimator.
    
    msesum one common MSE-optimal bandwidth selector for the sum of regression estimates (as opposed to difference thereof).
    
    msecomb1 for min(mserd,msesum).
    
    msecomb2 for median(msetwo,mserd,msesum), for each side of the cutoff separately.
    
    cerrd one common CER-optimal bandwidth selector for the RD treatment effect estimator.
    
    certwo two different CER-optimal bandwidth selectors (below and above the cutoff) for the RD treatment effect estimator.
    
    cersum one common CER-optimal bandwidth selector for the sum of regression estimates (as opposed to difference thereof).
    
    cercomb1 for min(cerrd,cersum).
    
    cercomb2 for median(certwo,cerrd,cersum), for each side of the cutoff separately.
    
    Note: MSE = Mean Square Error; CER = Coverage Error Rate. Default is bwselect=mserd. For details on implementation see Calonico, Cattaneo and Titiunik (2014a), Calonico, Cattaneo and Farrell (2018), and Calonico, Cattaneo, Farrell and Titiunik (2019), and the companion software articles.
    
    vce
    specifies the procedure used to compute the variance-covariance matrix estimator. Options are:

    nn for heteroskedasticity-robust nearest neighbor variance estimator with nnmatch the (minimum) number of neighbors to be used.

    hc0 for heteroskedasticity-robust plug-in residuals variance estimator without weights.

    hc1 for heteroskedasticity-robust plug-in residuals variance estimator with hc1 weights.

    hc2 for heteroskedasticity-robust plug-in residuals variance estimator with hc2 weights.

    hc3 for heteroskedasticity-robust plug-in residuals variance estimator with hc3 weights.

    cr1 for cluster-robust plug-in residuals variance estimator with degrees-of-freedom weights. Requires cluster.

    cr2 for the Bell-McCaffrey (2002) bias-reduced cluster-robust variance estimator (CRV2). Requires cluster.

    cr3 for the Pustejovsky-Tipton (2018) cluster-robust variance estimator (CRV3), approximately unbiased with small numbers of clusters. Requires cluster.

    The CR2/CR3 leverage correction applies to both the conventional and the robust bias-corrected standard errors, including when the point-estimation bandwidth h differs from the bias-correction bandwidth b; in that case the cluster leverage is computed from the bias (b) regression.

    Default is vce=nn without cluster, and vce=cr1 with cluster. When cluster is supplied and vce is a heteroskedasticity-robust option (nn/hc0/hc1/hc2/hc3), the estimator is automatically switched to the corresponding cluster-robust variant (cr1, cr2, or cr3) with a warning. If vce is a cluster option (cr1/cr2/cr3) but cluster is NULL, the estimator falls back to the corresponding heteroskedasticity-robust variant (hc1/hc2/hc3).

    cluster
    indicates the cluster ID variable used for cluster-robust variance estimation. When specified, only cluster-robust variance estimators (cr1, cr2, cr3) are valid; see the vce option for auto-switching rules.
    
    nnmatch	
    to be combined with for vce=nn for heteroskedasticity-robust nearest neighbor variance estimator with nnmatch indicating the minimum number of neighbors to be used. Default is nnmatch=3
    
    level	
    sets the confidence level for confidence intervals; default is level = 95.
    
    scalepar	
    specifies scaling factor for RD parameter of interest. This option is useful when the population parameter of interest involves a known multiplicative factor (e.g., sharp kink RD). Default is scalepar = 1 (no scaling).
    
    scaleregul	
    specifies scaling factor for the regularization term added to the denominator of the bandwidth selectors. Setting scaleregul = 0 removes the regularization term from the bandwidth selectors; default is scaleregul = 1.
    
    sharpbw	
    option to perform fuzzy RD estimation using a bandwidth selection procedure for the sharp RD model. This option is automatically selected if there is perfect compliance at either side of the cutoff.
    
    all	
    if specified, rdrobust reports three different procedures:
    
    (i) conventional RD estimates with conventional standard errors.
    
    (ii) bias-corrected estimates with conventional standard errors.
    
    (iii) bias-corrected estimates with robust standard errors.
    
    subset	
    an optional vector specifying a subset of observations to be used.
    
    masspoints	
    checks and controls for repeated observations in the running variable. Options are:
    
    (i) off: ignores the presence of mass points;
    
    (ii) check: looks for and reports the number of unique observations at each side of the cutoff.
    
    (iii) adjust: controls that the preliminary bandwidths used in the calculations contain a minimal number of unique observations. By default it uses 10 observations, but it can be manually adjusted with the option bwcheck.
    
    Default option is masspoints=adjust.
    
    bwcheck	
    if a positive integer is provided, the preliminary bandwidth used in the calculations is enlarged so that at least bwcheck unique observations are used.
    
    bwrestrict	
    if TRUE, computed bandwidths are restricted to lie within the range of x; default is bwrestrict = TRUE.
    
    stdvars
    if TRUE, x and y are standardized before computing the bandwidths; default is stdvars = FALSE.

    data
    optional pandas DataFrame. When supplied, `y`, `x`, `covs`, `cluster`,
    `fuzzy`, `weights`, and `subset` may be passed as column-name strings
    (or a list of names for `covs`), and a one-sided formula is accepted
    for `covs`.

    Returns
    -------
    N	
    vector with the sample sizes used to the left and to the right of the cutoff.
    
    N_h	
    vector with the effective sample sizes used to the left and to the right of the cutoff.
    
    c	
    cutoff value.
    
    p	
    order of the polynomial used for estimation of the regression function.
    
    q	
    order of the polynomial used for estimation of the bias of the regression function.
    
    bws	
    matrix containing the bandwidths used.
    
    tau_cl	
    conventional local-polynomial estimate to the left and to the right of the cutoff.
    
    tau_bc	
    bias-corrected local-polynomial estimate to the left and to the right of the cutoff.
    
    coef	
    vector containing conventional and bias-corrected local-polynomial RD estimates.
    
    se	
    vector containing conventional and robust standard errors of the local-polynomial RD estimates.
    
    bias	
    estimated bias for the local-polynomial RD estimator below and above the cutoff.
    
    beta_p_l	
    conventional p-order local-polynomial estimates to the left of the cutoff.
    
    beta_p_r	
    conventional p-order local-polynomial estimates to the right of the cutoff.
    
    V_cl_l	
    conventional variance-covariance matrix estimated below the cutoff.
    
    V_cl_r	
    conventional variance-covariance matrix estimated above the cutoff.
    
    V_rb_l	
    robust variance-covariance matrix estimated below the cutoff.
    
    V_rb_r	
    robust variance-covariance matrix estimated above the cutoff.
    
    z
    vector containing the z-statistics associated with conventional, bias-corrected and robust local-polynomial RD estimates.

    pv
    vector containing the p-values associated with conventional, bias-corrected and robust local-polynomial RD estimates.

    ci
    matrix containing the confidence intervals associated with conventional, bias-corrected and robust local-polynomial RD estimates.

    coef_covs
    coefficients of the additional covariates, only returned when `covs` is specified.

    Estimate
    1x4 matrix with columns `tau.us`, `tau.bc`, `se.us`, `se.rb`.

    N_b, M
    effective sample sizes for bias estimation (`N_b`) and number of unique observations (`M`) on the left and right of the cutoff.

    n_clust
    total number of clusters when `cluster` is specified; `None` otherwise. Per-side counts are in `n_clust_l` and `n_clust_r`.

    kernel, vce, bwselect, level, masspoints, rdmodel
    settings used (echoed for replay / printing).

    tau_T, se_T, z_T, pv_T, ci_T
    first-stage (treatment) inference for fuzzy designs: point estimate,
    standard error, z-statistic, p-value, and confidence interval for the
    conventional, bias-corrected, and robust versions of the T-component
    effect at the cutoff. `None` for sharp designs.

    beta_T_p_l, beta_T_p_r
    conventional p-order local-polynomial estimates of the T-component to
    the left and right of the cutoff (fuzzy designs only; `None` otherwise).


    See Also
    --------
    rdbwselect, rdplot
    
    Example
    ------- 
    >>> x = numpy.random.uniform(low=-1, high=1, size=1000)
    >>> y = 5+3*x+2*(x>=0) + numpy.random.uniform(size = 1000)
    >>> rdrobust(y,x)
    '''

    loc = resolve_data_args(data, y=y, x=x, cluster=cluster,
                            fuzzy=fuzzy, weights=weights, subset=subset)
    y, x      = loc['y'], loc['x']
    cluster   = loc['cluster']
    fuzzy     = loc['fuzzy']
    weights   = loc['weights']
    subset    = loc['subset']
    covs = resolve_covs(covs, data)

    # Check for errors in the INPUT
    if p is None and deriv is not None: p = deriv + 1
    if p is None: p = 1
    elif not np.isscalar(p) or p not in range(21):
        raise Exception('Polynomial order p incorrectly specified.')
        
    if q is None: q = p + 1
    elif not np.isscalar(q) or q not in range(21) or q<=p:
        raise Exception('Polynomial order (for bias correction) q incorrectly specified; q must be strictly greater than p.')
    
    if deriv is None: deriv = 0
    elif deriv not in range(21) or deriv>p:
        raise Exception('Derivative order incorrectly specified')
    
    kernel   = kernel.lower()
    kernel_list = ['uni','uniform','tri','triangular','epa','epanechnikov','']
    if kernel not in kernel_list:   
        raise Exception("kernel incorrectly specified")
    
    bwselect = bwselect.lower()
    bw_list = ['mserd','msetwo','msesum','msecomb1','msecomb2','cerrd',
               'certwo','cersum','cercomb1','cercomb2','']
    if bwselect not in bw_list:
        raise Exception("bwselect incorrectly specified")
    elif bwselect=="cct" or bwselect=="ik" or bwselect=="cv":
        raise Exception("bwselect options IK, CCT and CV have been deprecated. Please see help for new options")  
    
    vce = vce.lower()
    vce_list = ['nn','hc0','hc1','hc2','hc3','cr1','cr2','cr3','']
    if vce not in vce_list:
        raise Exception("vce incorrectly specified")

    # vce / cluster validation
    import warnings
    if cluster is not None:
        if vce in ('nn',''):
            vce = 'cr1'
        elif vce in ('hc0','hc1'):
            warnings.warn("vce='" + vce + "' is not a cluster option. Switching to vce='cr1'.")
            vce = 'cr1'
        elif vce == 'hc2':
            warnings.warn("vce='hc2' is not a cluster option. Switching to vce='cr2'.")
            vce = 'cr2'
        elif vce == 'hc3':
            warnings.warn("vce='hc3' is not a cluster option. Switching to vce='cr3'.")
            vce = 'cr3'
        elif vce not in ('cr1','cr2','cr3'):
            warnings.warn("vce='" + vce + "' is not a valid cluster option. Switching to vce='cr1'.")
            vce = 'cr1'
    else:
        if vce == 'cr1':
            warnings.warn("vce='cr1' requires cluster. Switching to vce='hc1'.")
            vce = 'hc1'
        elif vce == 'cr2':
            warnings.warn("vce='cr2' requires cluster. Switching to vce='hc2'.")
            vce = 'hc2'
        elif vce == 'cr3':
            warnings.warn("vce='cr3' requires cluster. Switching to vce='hc3'.")
            vce = 'hc3'
    
    if not (np.isscalar(level) and np.isfinite(level) and 0 < level < 100):
        raise Exception("level should be a single number in (0, 100)")

    if rho is not None and not (np.isscalar(rho) and np.isfinite(rho) and rho > 0):
        raise Exception("rho should be a single positive number")

    if bwcheck is not None:
        if not (np.isscalar(bwcheck) and np.isfinite(bwcheck) and bwcheck >= 1
                and bwcheck == round(bwcheck)):
            raise Exception("bwcheck must be a single positive integer")

    if (masspoints is not None and masspoints is not False
            and masspoints not in ("check", "adjust", "off", "")):
        raise Exception("masspoints must be one of 'check', 'adjust', 'off', or False")

    #=========================================================================
    # Tidy the Input and remove NAN

    x = np.array(x).reshape(-1,1)
    y = np.array(y).reshape(-1,1)

    # Validate auxiliary-vector lengths against length(x) BEFORE subset filtering,
    # so wrong-length inputs error explicitly instead of being silently recycled
    # by [subset] / [na_ok] indexing.
    n_orig = len(x)
    if len(y) != n_orig:
        raise Exception(f"'y' and 'x' must have equal length (got y={len(y)}, x={n_orig}).")
    if weights is not None and len(np.asarray(weights).reshape(-1)) != n_orig:
        raise Exception(f"'weights' must have length equal to length(x) (got {len(np.asarray(weights).reshape(-1))}, expected {n_orig}).")
    if fuzzy is not None and len(np.asarray(fuzzy).reshape(-1)) != n_orig:
        raise Exception(f"'fuzzy' must have length equal to length(x) (got {len(np.asarray(fuzzy).reshape(-1))}, expected {n_orig}).")
    if covs is not None:
        _covs_arr = np.asarray(covs)
        nc = _covs_arr.shape[0] if _covs_arr.ndim >= 1 else 0
        if nc != n_orig:
            raise Exception(f"'covs' must have nrow equal to length(x) (got {nc}, expected {n_orig}).")
    if cluster is not None and len(np.asarray(cluster).reshape(-1)) != n_orig:
        raise Exception(f"'cluster' must have length equal to length(x) (got {len(np.asarray(cluster).reshape(-1))}, expected {n_orig}).")
    if subset is not None:
        _subset_arr = np.asarray(subset)
        if _subset_arr.dtype == bool:
            if len(_subset_arr) != n_orig:
                raise Exception(f"Boolean 'subset' must have length equal to length(x) (got {len(_subset_arr)}, expected {n_orig}).")
        elif np.issubdtype(_subset_arr.dtype, np.integer) or np.issubdtype(_subset_arr.dtype, np.floating):
            if (not np.all(np.isfinite(_subset_arr)) or np.any(_subset_arr < 0)
                    or np.any(_subset_arr >= n_orig)
                    or not np.all(_subset_arr == _subset_arr.astype(int))):
                raise Exception(f"Numeric 'subset' must contain integer indices in 0..{n_orig - 1}.")
        else:
            raise Exception("'subset' must be boolean or integer.")

    if subset is not None:
        x = x[subset]
        y = y[subset]
        if len(x) == 0:
            raise Exception("'subset' removed all observations.")
    na_ok = complete_cases(x) & complete_cases(y)
    
    if cluster is not None:
        cluster = np.array(cluster).reshape(-1,1)
        if subset is not None:
            cluster = cluster[subset]
        if np.issubdtype(cluster.dtype, np.number):
            na_ok = na_ok & complete_cases(cluster)
    
    if covs is not None:
        try: covs_names = list(covs.columns)
        except: covs_names = ['z' + str(i+1) for i in range(ncol(covs))]         
        covs = np.array(covs).reshape(len(covs),-1)  
        covs_order = np.argsort([len(i) for i in covs_names])
        covs = covs[:,covs_order]
        covs_names = [covs_names[i] for i in covs_order]
        if subset is not None:
            covs = covs[subset,:]
        na_ok = na_ok & complete_cases(covs) 
    
    if fuzzy is not None:
        fuzzy = np.array(fuzzy).reshape(-1,1)
        if subset is not None:
            fuzzy = fuzzy[subset]
        na_ok = na_ok & complete_cases(fuzzy)
  
    if weights is not None:
        weights = np.array(weights).reshape(-1,1)
        if subset is not None:
            weights = weights[subset]
        na_ok = na_ok & complete_cases(weights) & (weights>=0).reshape(-1,)
    
    x = x[na_ok]
    y = y[na_ok]
    
    if covs is not None: covs = covs[na_ok,:]
    if fuzzy is not None: fuzzy   = fuzzy[na_ok]
    if cluster is not None: cluster = cluster[na_ok]
    if weights is not None: weights = weights[na_ok]
  
    if masspoints is None: masspoints = False

    if vce == "nn" or masspoints == "check" or masspoints == "adjust" : 
        order_x = np.argsort(x[:,0])
        x = x[order_x]
        y = y[order_x]
        if covs is not None: covs = covs[order_x,:]
        if fuzzy is not None: fuzzy   = fuzzy[order_x]
        if cluster is not None: cluster = cluster[order_x]
        if weights is not None: weights = weights[order_x]
    
    if c is None: c = 0
    if not (np.isscalar(c) and np.isreal(c) and np.isfinite(c)):
        raise Exception(f"Cutoff 'c' must be a single finite numeric value (received: {c}).")
    less_c = x<c
    X_l = x[less_c]
    X_r = x[np.invert(less_c)]
    Y_l = y[less_c]
    Y_r = y[np.invert(less_c)]
    x_min = np.min(x)  
    x_max = np.max(x)
    if c<=x_min or c>=x_max:
        raise Exception("c should be set within the range of x")
    range_l = np.abs(np.max(X_l)-np.min(X_l))
    range_r = np.abs(np.max(X_r)-np.min(X_r))
    N_l = len(X_l)
    N_r = len(X_r)
    N = N_r + N_l

    # Reject fully degenerate first stage (no variation, no jump). One-sided
    # non-compliance falls through to the perf_comp branch downstream.
    if fuzzy is not None:
        T_l_chk = fuzzy[less_c]
        T_r_chk = fuzzy[np.invert(less_c)]
        if (np.var(T_l_chk) == 0 and np.var(T_r_chk) == 0
                and abs(np.mean(T_l_chk) - np.mean(T_r_chk)) < np.sqrt(np.finfo(float).eps)):
            raise ValueError(
                "Fuzzy RD: first-stage variable has no variation and no jump "
                "at the cutoff. The fuzzy estimator is not identified."
            )
    quant = -sct.norm.ppf(q = np.abs((1-(level/100))/2))
    
    # Internal mapping: cr1/cr2/cr3 -> internal vce string + flags
    vce_type = "NN"
    if vce=="cr1": vce_type = "CR1"; vce = "hc1"
    elif vce=="cr2": vce_type = "CR2"; vce = "crv2"
    elif vce=="cr3": vce_type = "CR3"; vce = "crv3"
    elif vce=="hc0": vce_type = "HC0"
    elif vce=="hc1": vce_type = "HC1"
    elif vce=="hc2": vce_type = "HC2"
    elif vce=="hc3": vce_type = "HC3"
    crv3 = (vce=="crv3") and cluster is not None
    crv2 = (vce=="crv2") and cluster is not None

    #===================================================================================
    # Check for COLLINEARITY
    
    covs_drop_coll = dZ = 0
    if covs_drop: covs_drop_coll = 1 
    if covs is not None:
        dZ = ncol(covs)
        if covs_drop: 
            covs_check = covs_drop_fun(covs)
            if ncol(covs_check) < dZ:
                covs = covs_check
                dZ = ncol(covs)
    # ================================================================
                  
    if h is not None: 
        bwselect = "Manual"
        if rho is None:
            rho = 1
            if b is None:
                b = h
        else:
            b = h/rho
            
    if N<20:
        print("Not enough observations to perform bandwidth calculations. Estimates computed using entire sample")
        h = b = np.max(range_l,range_r)
        bwselect = "Manual"
  
    if kernel=="epanechnikov" or kernel=="epa":
        kernel_type = "Epanechnikov"
    elif kernel=="uniform" or kernel=="uni":
        kernel_type = "Uniform"
    else:
        kernel_type = "Triangular"

    # vce_type already set above; no need to re-derive

    M_l = N_l
    M_r = N_r
    X_uniq_l = X_uniq_r = None
    # X_uniq_l/r and M_l/M_r are needed by both the masspoints inspection
    # branch (mass detection) AND the bwcheck branch (bw_min_l/r). Hoist
    # the unique-x computation so bwcheck works even when masspoints='off'.
    if masspoints=="check" or masspoints=="adjust" or bwcheck is not None:
      X_uniq_l = np.sort(np.unique(X_l))[::-1]
      X_uniq_r = np.unique(X_r)
      M_l = len(X_uniq_l)
      M_r = len(X_uniq_r)
    if masspoints=="check" or masspoints=="adjust":
      mass_l = 1-M_l/N_l
      mass_r = 1-M_r/N_r
      if mass_l>=0.2 or mass_r>=0.2:
        print("Mass points detected in the running variable.")
        if masspoints=="check": print("Try using option masspoints=adjust.")
        if bwcheck is None and masspoints=="adjust": bwcheck = 10

    ######### Calculate bandwidth

    if h is None:
        # Inline bandwidth computation: skip rdbwselect() entrypoint (which
        # redoes NA/sort/split/dups prep we already have) and call the shared
        # compute helper directly with our already-prepared data.

        # Rescaling (mirrors rdbwselect's top block)
        x_q25_bw, x_q75_bw = quantile_type2(x, [0.25, 0.75])
        x_iq_bw = x_q75_bw - x_q25_bw
        BWp = min(np.std(x, ddof=1), x_iq_bw / 1.349)
        x_sd_bw = 1.0
        y_sd_bw = 1.0
        c_bw_arg = c
        X_l_bw, X_r_bw = X_l, X_r
        Y_l_bw, Y_r_bw = Y_l, Y_r
        x_min_bw, x_max_bw = x_min, x_max
        # Range for BW uses distance from cutoff to the extreme on each side
        # (NOT the running-variable-width that rdrobust caches for its own use).
        range_l_bw = float(np.abs(c - np.min(X_l)))
        range_r_bw = float(np.abs(c - np.max(X_r)))
        if stdvars:
            y_sd_bw = float(np.std(y, ddof=1))
            x_sd_bw = float(np.std(x, ddof=1))
            X_l_bw = X_l / x_sd_bw
            X_r_bw = X_r / x_sd_bw
            Y_l_bw = Y_l / y_sd_bw
            Y_r_bw = Y_r / y_sd_bw
            c_bw_arg = c / x_sd_bw
            x_min_bw = x_min / x_sd_bw
            x_max_bw = x_max / x_sd_bw
            range_l_bw = range_l_bw / x_sd_bw
            range_r_bw = range_r_bw / x_sd_bw
            BWp = min(1.0, (x_iq_bw / x_sd_bw) / 1.349)

        # Kernel -> C_c
        if kernel == "epanechnikov" or kernel == "epa":
            C_c_bw = 2.34
        elif kernel == "uniform" or kernel == "uni":
            C_c_bw = 1.843
        else:
            C_c_bw = 2.576

        # Full-side data for BW (not ind-filtered).
        T_l_bw = T_r_bw = None
        if fuzzy is not None:
            T_l_bw_raw = fuzzy[x < c]
            T_r_bw_raw = fuzzy[x >= c]
            if np.var(T_l_bw_raw) == 0 or np.var(T_r_bw_raw) == 0 or sharpbw:
                pass  # T stays None for sharp-style BW
            else:
                T_l_bw, T_r_bw = T_l_bw_raw, T_r_bw_raw
        Z_l_bw = Z_r_bw = None
        if covs is not None:
            Z_l_bw = covs[(x < c).reshape(-1), :]
            Z_r_bw = covs[(x >= c).reshape(-1), :]
        C_l_bw = C_r_bw = None
        g_l_bw = g_r_bw = None
        if cluster is not None:
            C_l_bw = cluster[x < c]
            C_r_bw = cluster[x >= c]
            g_l_bw = len(np.unique(C_l_bw))
            g_r_bw = len(np.unique(C_r_bw))
        fw_l_bw = fw_r_bw = 0
        if weights is not None:
            fw_l_bw = weights[x < c]
            fw_r_bw = weights[x >= c]

        # dups for nn (full-sample, not effective-sample)
        dups_l_bw = np.zeros(N_l, dtype=int)
        dupsid_l_bw = np.zeros(N_l, dtype=int)
        dups_r_bw = np.zeros(N_r, dtype=int)
        dupsid_r_bw = np.zeros(N_r, dtype=int)
        if vce == "nn":
            aux_l = pd.DataFrame({'nn_l': np.ones(N_l), 'X_l': X_l_bw})
            dups_l_bw = aux_l.groupby('X_l')['nn_l'].transform('sum').values.astype(int)
            dupsid_l_bw = aux_l.groupby('X_l')['nn_l'].transform('cumsum').values.astype(int)
            aux_r = pd.DataFrame({'nn_r': np.ones(N_r), 'X_r': X_r_bw})
            dups_r_bw = aux_r.groupby('X_r')['nn_r'].transform('sum').values.astype(int)
            dupsid_r_bw = aux_r.groupby('X_r')['nn_r'].transform('cumsum').values.astype(int)

        # Apply internal vce mapping (cr1/crv2/crv3) to match rdbwselect
        vce_bw = vce
        if vce_bw == "cr1": vce_bw = "hc1"
        elif vce_bw == "cr2": vce_bw = "crv2"
        elif vce_bw == "cr3": vce_bw = "crv3"

        bws_df, _ = _rdbwselect_compute(
            Y_l=Y_l_bw, Y_r=Y_r_bw, X_l=X_l_bw, X_r=X_r_bw,
            T_l=T_l_bw, T_r=T_r_bw, Z_l=Z_l_bw, Z_r=Z_r_bw,
            C_l=C_l_bw, C_r=C_r_bw, fw_l=fw_l_bw, fw_r=fw_r_bw,
            dups_l=dups_l_bw, dups_r=dups_r_bw,
            dupsid_l=dupsid_l_bw, dupsid_r=dupsid_r_bw,
            N_l=N_l, N_r=N_r, N=N, M_l=M_l, M_r=M_r,
            M=(M_l + M_r) if masspoints in ("check", "adjust") else N,
            X_uniq_l=X_uniq_l, X_uniq_r=X_uniq_r,
            x_min=x_min_bw, x_max=x_max_bw,
            range_l=range_l_bw, range_r=range_r_bw, x_sd=x_sd_bw,
            c=c_bw_arg, p=p, q=q, deriv=deriv, kernel=kernel,
            C_c=C_c_bw, BWp=BWp,
            bwselect=bwselect, scaleregul=scaleregul,
            vce=vce_bw, nnmatch=nnmatch, covs_drop_coll=covs_drop_coll,
            bwcheck=bwcheck, bwrestrict=bwrestrict, masspoints=masspoints,
            g_l=g_l_bw, g_r=g_r_bw, cluster_present=(cluster is not None),
            all=False,
        )
        h_l, h_r, b_l, b_r = bws_df.iloc[0].values
        if rho is not None:
            b_l = h_l / rho
            b_r = h_r / rho
    else:
        if np.isscalar(h): h_l = h_r = h
        elif len(h)==2: h_l, h_r = h
        if b is None:      
            b_l = h_l
            b_r = h_r
        else:
          if np.isscalar(b): b_l = b_r = b
          elif len(b)==2: b_l, b_r = b
    w_h_l = rdrobust_kweight(X_l,c,h_l,kernel)	
    w_h_r = rdrobust_kweight(X_r,c,h_r,kernel)
    w_b_l = rdrobust_kweight(X_l,c,b_l,kernel)	
    w_b_r = rdrobust_kweight(X_r,c,b_r,kernel)
    
    if weights is not None:
        fw_l = weights[x<c]  
        fw_r = weights[x>=c]
        w_h_l = fw_l*w_h_l
        w_h_r = fw_r*w_h_r
        w_b_l = fw_l*w_b_l
        w_b_r = fw_r*w_b_r			
      
    ind_h_l = w_h_l> 0
    ind_h_r = w_h_r> 0
    ind_b_l = w_b_l> 0
    ind_b_r = w_b_r> 0
    N_h_l = np.sum(ind_h_l)
    N_b_l = np.sum(ind_b_l)
    N_h_r = np.sum(ind_h_r)
    N_b_r = np.sum(ind_b_r)
    
    ind_l = ind_b_l 
    ind_r = ind_b_r
    if h_l>b_l: ind_l = ind_h_l   
    if h_r>b_r: ind_r = ind_h_r   
    
    eN_l = np.sum(ind_l)
    eN_r = np.sum(ind_r)
    eY_l  = Y_l[ind_l]
    eY_r  = Y_r[ind_r]
    eX_l  = X_l[ind_l]
    eX_r  = X_r[ind_r]
    W_h_l = w_h_l[ind_l].reshape(-1,1)
    W_h_r = w_h_r[ind_r].reshape(-1,1)
    W_b_l = w_b_l[ind_l].reshape(-1,1)
    W_b_r = w_b_r[ind_r].reshape(-1,1)
    
    edups_l = np.zeros(eN_l).astype(int)
    edupsid_l = np.zeros(eN_l).astype(int)
    edups_r = np.zeros(eN_r).astype(int)
    edupsid_r = np.zeros(eN_r).astype(int)
    if vce=="nn":
        aux_l  = pd.DataFrame({'nn_l': np.ones(eN_l), 'eX_l': eX_l })
        edups_l   = aux_l.groupby('eX_l')['nn_l'].transform('sum').values.astype(int)
        edupsid_l = aux_l.groupby('eX_l')['nn_l'].transform('cumsum').values.astype(int)

        aux_r  = pd.DataFrame({'nn_r': np.ones(eN_r), 'eX_r': eX_r })
        edups_r   = aux_r.groupby('eX_r')['nn_r'].transform('sum').values.astype(int)
        edupsid_r = aux_r.groupby('eX_r')['nn_r'].transform('cumsum').values.astype(int)

    u_l = ((eX_l-c)/h_l).reshape(-1,1)
    u_r = ((eX_r-c)/h_r).reshape(-1,1)
    # Q1: Vandermonde via successive multiplication (built in funs._vander).
    from .funs import _vander
    R_q_l = _vander(eX_l - c, q)
    R_q_r = _vander(eX_r - c, q)
    R_p_l = R_q_l[:,:(p+1)]
    R_p_r = R_q_r[:,:(p+1)]

    # Precompute kernel-weighted design matrices (used 6+ times each below)
    RW_p_l = R_p_l * W_h_l
    RW_p_r = R_p_r * W_h_r
    RW_q_l_b = R_q_l * W_b_l
    RW_q_r_b = R_q_r * W_b_r

    # Computing RD estimates

    L_l = crossprod(RW_p_l,u_l**(p+1))
    L_r = crossprod(RW_p_r,u_r**(p+1))
    # Py-2: cache G alongside invG so CRV3 branches can skip the inv(invG) round-trip.
    sqW_q_l_R = np.sqrt(W_b_l).reshape(-1,1) * R_q_l
    sqW_q_r_R = np.sqrt(W_b_r).reshape(-1,1) * R_q_r
    G_q_l     = crossprod(sqW_q_l_R)
    G_q_r     = crossprod(sqW_q_r_R)
    invG_q_l  = inv_chol(G_q_l)
    invG_q_r  = inv_chol(G_q_r)
    invG_p_l  = qrXXinv((np.sqrt(W_h_l)*R_p_l))
    invG_p_r  = qrXXinv((np.sqrt(W_h_r)*R_p_r))
    e_p1 = np.zeros((q+1,1))
    e_p1[p+1] = 1
    e_v  = np.zeros((p+1,1))
    e_v[deriv] = 1
    # Q_q mathematically reduces to: RW_p - h^(p+1) * outer(m, L), where
    # m[i] = (R_q @ invG_q[p+1, :]) * W_b (row p+1 is the only one that matters
    # because e_p1 selects it). This avoids a (q+1, n) intermediate and a
    # zero-col-heavy (p+1, q+1) matmul.
    m_l = (R_q_l @ invG_q_l[p+1, :]).reshape(-1, 1) * W_b_l
    m_r = (R_q_r @ invG_q_r[p+1, :]).reshape(-1, 1) * W_b_r
    Q_q_l = RW_p_l - h_l**(p+1) * m_l * L_l.reshape(1, -1)
    Q_q_r = RW_p_r - h_r**(p+1) * m_r * L_r.reshape(1, -1)
    D_l = eY_l.copy()
    D_r = eY_r.copy()
    eC_l = eC_r = eT_l = eT_r = eZ_l = eZ_r = None
    
    dT = 0
    if fuzzy is not None:
        dT = 1
        T_l  = fuzzy[x<c]
        eT_l  = T_l[ind_l]
        T_r  = fuzzy[x>=c]
        eT_r  = T_r[ind_r]
        D_l  = np.column_stack((D_l,eT_l))
        D_r = np.column_stack((D_r,eT_r))
  
    if covs is not None:
        Z_l  = covs[(x<c).reshape(-1),:]
        eZ_l = Z_l[ind_l,:]
        Z_r  = covs[(x>=c).reshape(-1),:]
        eZ_r = Z_r[ind_r,:]
        D_l  = np.column_stack((D_l,eZ_l))
        D_r = np.column_stack((D_r,eZ_r))
        U_p_l = crossprod(RW_p_l,D_l)
        U_p_r = crossprod(RW_p_r,D_r)
  
    cidx_l = cidx_r = None
    if cluster is not None:
        C_l  = cluster[x<c]
        C_r = cluster[x>=c]
        eC_l  = C_l[ind_l]
        eC_r  = C_r[ind_r]
        # T3: precompute per-side cluster index once for all V_Y/V_T vce calls.
        from .funs import _build_cluster_idx
        cidx_l = _build_cluster_idx(np.asarray(eC_l).reshape(-1))
        cidx_r = _build_cluster_idx(np.asarray(eC_r).reshape(-1))

    
    beta_p_l = tomat(np.matmul(invG_p_l,crossprod(RW_p_l,D_l)))
    beta_q_l = tomat(np.matmul(invG_q_l,crossprod(RW_q_l_b,D_l)))
    beta_bc_l = tomat(np.matmul(invG_p_l,crossprod(Q_q_l,D_l)))
    beta_p_r = tomat(np.matmul(invG_p_r,crossprod(RW_p_r,D_r)))
    beta_q_r = tomat(np.matmul(invG_q_r,crossprod(RW_q_r_b,D_r)))
    beta_bc_r = tomat(np.matmul(invG_p_r,crossprod(Q_q_r,D_r)))
    beta_p  = beta_p_r  - beta_p_l
    beta_bc = beta_bc_r - beta_bc_l
    

    if covs is None:
        tau_cl = scalepar*math.factorial(deriv)*beta_p[deriv,0]
        tau_Y_cl = tau_cl.copy()
        tau_bc = scalepar*math.factorial(deriv)*beta_bc[deriv,0]
        tau_Y_bc = tau_bc.copy()
        s_Y = 1        
        tau_Y_cl_l = scalepar*math.factorial(deriv)*beta_p_l[deriv,0]
        tau_Y_cl_r = scalepar*math.factorial(deriv)*beta_p_r[deriv,0]
        tau_Y_bc_l = scalepar*math.factorial(deriv)*beta_bc_l[deriv,0]
        tau_Y_bc_r = scalepar*math.factorial(deriv)*beta_bc_r[deriv,0]
        bias_l = tau_Y_cl_l-tau_Y_bc_l
        bias_r = tau_Y_cl_r-tau_Y_bc_r 
      
        if fuzzy is not None:
            tau_T_cl = math.factorial(deriv)*beta_p[deriv,1]
            tau_T_bc = math.factorial(deriv)*beta_bc[deriv,1]
            tau_cl   = tau_Y_cl/tau_T_cl
            s_Y      = np.array([1/tau_T_cl , -(tau_Y_cl/tau_T_cl**2)])
            B_F      = np.array([tau_Y_cl-tau_Y_bc , tau_T_cl-tau_T_bc])
            tau_bc   = tau_cl - np.matmul(s_Y.T,B_F)
            sV_T     = np.array([0, 1])

            tau_T_cl_l = math.factorial(deriv)*beta_p_l[deriv,1]
            tau_T_cl_r = math.factorial(deriv)*beta_p_r[deriv,1]
            tau_T_bc_l = math.factorial(deriv)*beta_bc_l[deriv,1]
            tau_T_bc_r = math.factorial(deriv)*beta_bc_r[deriv,1]
            B_F_l = np.array([tau_Y_cl_l-tau_Y_bc_l, tau_T_cl_l-tau_T_bc_l])
            B_F_r = np.array([tau_Y_cl_r-tau_Y_bc_r, tau_T_cl_r-tau_T_bc_r])
            bias_l = np.matmul(s_Y.T,B_F_l)
            bias_r = np.matmul(s_Y.T,B_F_r)

            beta_T_p_l = scalepar*math.factorial(deriv)*beta_p_l[:,1]
            beta_T_p_r = scalepar*math.factorial(deriv)*beta_p_r[:,1]
    else:
        ZWD_p_l  = crossprod(eZ_l*W_h_l,D_l)
        ZWD_p_r  = crossprod(eZ_r*W_h_r,D_r)
        colsZ = np.arange(1+dT, max(2+dT+dZ-1,2+dT))
        UiGU_p_l =  crossprod(U_p_l[:,colsZ],np.matmul(invG_p_l,U_p_l)) 
        UiGU_p_r =  crossprod(U_p_r[:,colsZ],np.matmul(invG_p_r,U_p_r))  
        ZWZ_p_l = ZWD_p_l[:,colsZ] - UiGU_p_l[:,colsZ] 
        ZWZ_p_r = ZWD_p_r[:,colsZ] - UiGU_p_r[:,colsZ]     
        ZWY_p_l = ZWD_p_l[:,:(1+dT)] - UiGU_p_l[:,:(1+dT)] 
        ZWY_p_r = ZWD_p_r[:,:(1+dT)] - UiGU_p_r[:,:(1+dT)]      
        ZWZ_p = ZWZ_p_r + ZWZ_p_l
        ZWY_p = ZWY_p_r + ZWY_p_l
        if covs_drop_coll == 0:
            gamma_p = np.matmul(inv_chol(ZWZ_p),ZWY_p)
        elif covs_drop_coll == 1:
            gamma_p = np.matmul(np.linalg.pinv(ZWZ_p),ZWY_p)
        s_Y = (np.hstack([1,-gamma_p[:,0]])).reshape(-1,1)
        
        if fuzzy is None:
            tau_cl = np.matmul(scalepar*s_Y.T,beta_p[deriv,:]).item()
            tau_bc = np.matmul(scalepar*s_Y.T,beta_bc[deriv,:]).item()
            tau_Y_cl_l = np.matmul(scalepar*s_Y.T,beta_p_l[deriv,:]).item()
            tau_Y_cl_r = np.matmul(scalepar*s_Y.T,beta_p_r[deriv,:]).item()
            tau_Y_bc_l = np.matmul(scalepar*s_Y.T,beta_bc_l[deriv,:]).item()
            tau_Y_bc_r = np.matmul(scalepar*s_Y.T,beta_bc_r[deriv,:]).item()
            bias_l = tau_Y_cl_l - tau_Y_bc_l
            bias_r = tau_Y_cl_r - tau_Y_bc_r
        else:
            s_T  = np.append(1,-gamma_p[:,1])
            sV_T = np.concatenate(([0, 1], -gamma_p[:,1]))

            tau_Y_cl = np.matmul(scalepar*math.factorial(deriv)*s_Y.T,
                        np.append(beta_p[deriv,0], beta_p[deriv,colsZ]).reshape(-1,1)).item()

            tau_T_cl = np.matmul(math.factorial(deriv)*s_T.T,
                        np.append(beta_p[deriv,1], beta_p[deriv,colsZ]).reshape(-1,1)).item()

            tau_Y_bc = np.matmul(scalepar*math.factorial(deriv)*s_Y.T,
                        np.append(beta_bc[deriv,0], beta_bc[deriv,colsZ]).reshape(-1,1)).item()

            tau_T_bc = np.matmul(math.factorial(deriv)*s_T.T,
                        np.append(beta_bc[deriv,1], beta_bc[deriv,colsZ]).reshape(-1,1)).item()

            tau_Y_cl_l = np.matmul(scalepar*math.factorial(deriv)*s_Y.T,
                        np.append(beta_p_l[deriv,0], beta_p_l[deriv,colsZ]).reshape(-1,1)).item()

            tau_Y_cl_r = np.matmul(scalepar*math.factorial(deriv)*s_Y.T,
                        np.append(beta_p_r[deriv,0], beta_p_r[deriv,colsZ]).reshape(-1,1)).item()

            tau_Y_bc_l = np.matmul(scalepar*math.factorial(deriv)*s_Y.T,
                        np.append(beta_bc_l[deriv,0], beta_bc_l[deriv,colsZ]).reshape(-1,1)).item()

            tau_Y_bc_r = np.matmul(scalepar*math.factorial(deriv)*s_Y.T,
                        np.append(beta_bc_r[deriv,0], beta_bc_r[deriv,colsZ]).reshape(-1,1)).item()

            tau_T_cl_l = np.matmul(math.factorial(deriv)*s_T.T,
                        np.append(beta_p_l[deriv,1], beta_p_l[deriv,colsZ]).reshape(-1,1)).item()

            tau_T_cl_r = np.matmul(math.factorial(deriv)*s_T.T,
                        np.append(beta_p_r[deriv,1], beta_p_r[deriv,colsZ]).reshape(-1,1)).item()

            tau_T_bc_l = np.matmul(math.factorial(deriv)*s_T.T,
                        np.append(beta_bc_l[deriv,1], beta_bc_l[deriv,colsZ]).reshape(-1,1)).item()

            tau_T_bc_r = np.matmul(math.factorial(deriv)*s_T.T,
                        np.append(beta_bc_r[deriv,1], beta_bc_r[deriv,colsZ]).reshape(-1,1)).item()

            beta_T_p_l = math.factorial(deriv) * (np.column_stack([beta_p_l[:,1], beta_p_l[:,colsZ]]) @ s_T)
            beta_T_p_r = math.factorial(deriv) * (np.column_stack([beta_p_r[:,1], beta_p_r[:,colsZ]]) @ s_T)

            tau_cl = tau_Y_cl/tau_T_cl
            B_F   = np.array([tau_Y_cl-tau_Y_bc, tau_T_cl-tau_T_bc])
            B_F_l = np.array([tau_Y_cl_l-tau_Y_bc_l, tau_T_cl_l-tau_T_bc_l])
            B_F_r = np.array([tau_Y_cl_r-tau_Y_bc_r, tau_T_cl_r-tau_T_bc_r])
            
            s_Y = np.array([1/tau_T_cl , -(tau_Y_cl/tau_T_cl**2)])
            tau_bc = tau_cl - np.dot(s_Y,B_F)
            
            bias_l = np.dot(s_Y,B_F_l)
            bias_r = np.dot(s_Y,B_F_r)
          
            s_Y = np.hstack([s_Y,-(1/tau_T_cl)*gamma_p[:,0] + (tau_Y_cl/tau_T_cl**2)*gamma_p[:,1]])
 
    ####### Computing variance-covariance matrix.

    hii_p_l = hii_p_r = hii_q_l = hii_q_r = 0
    predicts_p_l = predicts_p_r = predicts_q_l = predicts_q_r = 0
    if (vce in ("hc0","hc1","hc2","hc3")) and not crv3 and not crv2:
        predicts_p_l = np.matmul(R_p_l,beta_p_l)
        predicts_p_r = np.matmul(R_p_r,beta_p_r)
        predicts_q_l = np.matmul(R_q_l,beta_q_l)
        predicts_q_r = np.matmul(R_q_r,beta_q_r)
        if vce=="hc2" or vce=="hc3":
            hii_p_l = np.sum(np.matmul(R_p_l,invG_p_l)*RW_p_l, axis = 1)
            hii_p_r = np.sum(np.matmul(R_p_r,invG_p_r)*RW_p_r, axis = 1)
            hii_q_l = np.sum(np.matmul(R_q_l,invG_q_l)*RW_q_l_b, axis = 1)
            hii_q_r = np.sum(np.matmul(R_q_r,invG_q_r)*RW_q_r_b, axis = 1)
    elif crv3 or crv2:
        # CRV2/CRV3: compute predictions for raw residuals
        predicts_p_l = np.matmul(R_p_l,beta_p_l)
        predicts_p_r = np.matmul(R_p_r,beta_p_r)
        predicts_q_l = np.matmul(R_q_l,beta_q_l)
        predicts_q_r = np.matmul(R_q_r,beta_q_r)

    res_h_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_p_l, hii_p_l, vce, nnmatch, edups_l, edupsid_l, p+1, crv3=crv3, crv2=crv2, has_cluster=cluster is not None)
    res_h_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_p_r, hii_p_r, vce, nnmatch, edups_r, edupsid_r, p+1, crv3=crv3, crv2=crv2, has_cluster=cluster is not None)
    if vce=="nn":
        res_b_l = res_h_l
        res_b_r = res_h_r
    else:
        res_b_l = rdrobust_res(eX_l, eY_l, eT_l, eZ_l, predicts_q_l, hii_q_l, vce, nnmatch, edups_l, edupsid_l, q+1, crv3=crv3, crv2=crv2, has_cluster=cluster is not None)
        res_b_r = rdrobust_res(eX_r, eY_r, eT_r, eZ_r, predicts_q_r, hii_q_r, vce, nnmatch, edups_r, edupsid_r, q+1, crv3=crv3, crv2=crv2, has_cluster=cluster is not None)

    # CRV2/CRV3: sqrt-kernel-weighted design matrices for symmetric hat-matrix blocks
    sqrtRX_p_l = sqrtRX_p_r = None
    if crv3 or crv2:
        sqrtRX_p_l = R_p_l * np.sqrt(W_h_l)
        sqrtRX_p_r = R_p_r * np.sqrt(W_h_r)
    crv_invG_l = invG_p_l if (crv3 or crv2) else None
    crv_invG_r = invG_p_r if (crv3 or crv2) else None

    V_Y_cl_l = np.matmul(invG_p_l,np.matmul(rdrobust_vce(dT+dZ, s_Y, RW_p_l, res_h_l, eC_l, cluster_idx=cidx_l, invG=crv_invG_l, sqrtRX=sqrtRX_p_l, crv2=crv2),invG_p_l))
    V_Y_cl_r = np.matmul(invG_p_r,np.matmul(rdrobust_vce(dT+dZ, s_Y, RW_p_r, res_h_r, eC_r, cluster_idx=cidx_r, invG=crv_invG_r, sqrtRX=sqrtRX_p_r, crv2=crv2),invG_p_r))

    # V_rb: when h == b, use q-regression formulation directly (clean CRV2/CRV3).
    # When h != b with CRV2/CRV3, use Q_q sandwich + q-regression cluster leverage
    # (matrix analog of hii_q for HC2/HC3). Otherwise CR1.
    hb_match = (h_l == b_l) and (h_r == b_r)
    if hb_match and (crv3 or crv2):
        sqrtRX_q_l = R_q_l * np.sqrt(W_h_l)
        sqrtRX_q_r = R_q_r * np.sqrt(W_h_r)
        RW_q_l_h = R_q_l * W_h_l
        RW_q_r_h = R_q_r * W_h_r
        V_Y_rb_l = np.matmul(invG_q_l,np.matmul(rdrobust_vce(dT+dZ, s_Y, RW_q_l_h, res_b_l, eC_l, cluster_idx=cidx_l, invG=invG_q_l, sqrtRX=sqrtRX_q_l, crv2=crv2, G=(None if crv2 else G_q_l)),invG_q_l))
        V_Y_rb_r = np.matmul(invG_q_r,np.matmul(rdrobust_vce(dT+dZ, s_Y, RW_q_r_h, res_b_r, eC_r, cluster_idx=cidx_r, invG=invG_q_r, sqrtRX=sqrtRX_q_r, crv2=crv2, G=(None if crv2 else G_q_r)),invG_q_r))
    elif hb_match and cluster is not None:
        V_Y_rb_l = np.matmul(invG_q_l,np.matmul(rdrobust_vce(dT+dZ, s_Y, R_q_l*W_h_l, res_b_l, eC_l, cluster_idx=cidx_l),invG_q_l))
        V_Y_rb_r = np.matmul(invG_q_r,np.matmul(rdrobust_vce(dT+dZ, s_Y, R_q_r*W_h_r, res_b_r, eC_r, cluster_idx=cidx_r),invG_q_r))
    elif crv3 or crv2:
        M_rb_l = rdrobust_vce_qq_cluster(Q_q_l, R_q_l, W_b_l, invG_q_l, res_b_l, eC_l, cluster_idx=cidx_l, crv2=crv2, s=s_Y, d=dT+dZ, G_q=(None if crv2 else G_q_l))
        M_rb_r = rdrobust_vce_qq_cluster(Q_q_r, R_q_r, W_b_r, invG_q_r, res_b_r, eC_r, cluster_idx=cidx_r, crv2=crv2, s=s_Y, d=dT+dZ, G_q=(None if crv2 else G_q_r))
        V_Y_rb_l = invG_p_l @ M_rb_l @ invG_p_l
        V_Y_rb_r = invG_p_r @ M_rb_r @ invG_p_r
    else:
        # k_override=q+1 aligns the CR1 df correction with the q-regression
        # path used at h=b, so SE is continuous across the h=b boundary.
        # Non-cluster paths ignore k (df adjustment is in rdrobust_res's w factor).
        V_Y_rb_l = np.matmul(invG_p_l,np.matmul(rdrobust_vce(dT+dZ, s_Y, Q_q_l, res_b_l, eC_l, cluster_idx=cidx_l, k_override=q+1),invG_p_l))
        V_Y_rb_r = np.matmul(invG_p_r,np.matmul(rdrobust_vce(dT+dZ, s_Y, Q_q_r, res_b_r, eC_r, cluster_idx=cidx_r, k_override=q+1),invG_p_r))

    # First-stage (T) variance for fuzzy designs: uses sV_T selector (picks T coefficient
    # net of partialled-out Z), not the delta-method s_Y. Mirrors V_Y_rb branch structure.
    if fuzzy is not None:
        V_T_cl_l = np.matmul(invG_p_l,np.matmul(rdrobust_vce(dT+dZ, sV_T, RW_p_l, res_h_l, eC_l, cluster_idx=cidx_l, invG=crv_invG_l, sqrtRX=sqrtRX_p_l, crv2=crv2),invG_p_l))
        V_T_cl_r = np.matmul(invG_p_r,np.matmul(rdrobust_vce(dT+dZ, sV_T, RW_p_r, res_h_r, eC_r, cluster_idx=cidx_r, invG=crv_invG_r, sqrtRX=sqrtRX_p_r, crv2=crv2),invG_p_r))
        if hb_match and (crv3 or crv2):
            V_T_rb_l = np.matmul(invG_q_l,np.matmul(rdrobust_vce(dT+dZ, sV_T, RW_q_l_h, res_b_l, eC_l, cluster_idx=cidx_l, invG=invG_q_l, sqrtRX=sqrtRX_q_l, crv2=crv2, G=(None if crv2 else G_q_l)),invG_q_l))
            V_T_rb_r = np.matmul(invG_q_r,np.matmul(rdrobust_vce(dT+dZ, sV_T, RW_q_r_h, res_b_r, eC_r, cluster_idx=cidx_r, invG=invG_q_r, sqrtRX=sqrtRX_q_r, crv2=crv2, G=(None if crv2 else G_q_r)),invG_q_r))
        elif hb_match and cluster is not None:
            V_T_rb_l = np.matmul(invG_q_l,np.matmul(rdrobust_vce(dT+dZ, sV_T, R_q_l*W_h_l, res_b_l, eC_l, cluster_idx=cidx_l),invG_q_l))
            V_T_rb_r = np.matmul(invG_q_r,np.matmul(rdrobust_vce(dT+dZ, sV_T, R_q_r*W_h_r, res_b_r, eC_r, cluster_idx=cidx_r),invG_q_r))
        elif crv3 or crv2:
            M_T_rb_l = rdrobust_vce_qq_cluster(Q_q_l, R_q_l, W_b_l, invG_q_l, res_b_l, eC_l, cluster_idx=cidx_l, crv2=crv2, s=sV_T, d=dT+dZ, G_q=(None if crv2 else G_q_l))
            M_T_rb_r = rdrobust_vce_qq_cluster(Q_q_r, R_q_r, W_b_r, invG_q_r, res_b_r, eC_r, cluster_idx=cidx_r, crv2=crv2, s=sV_T, d=dT+dZ, G_q=(None if crv2 else G_q_r))
            V_T_rb_l = invG_p_l @ M_T_rb_l @ invG_p_l
            V_T_rb_r = invG_p_r @ M_T_rb_r @ invG_p_r
        else:
            # See V_Y_rb path-B comment above (k_override aligns CR1 df with h=b).
            V_T_rb_l = np.matmul(invG_p_l,np.matmul(rdrobust_vce(dT+dZ, sV_T, Q_q_l, res_b_l, eC_l, cluster_idx=cidx_l, k_override=q+1),invG_p_l))
            V_T_rb_r = np.matmul(invG_p_r,np.matmul(rdrobust_vce(dT+dZ, sV_T, Q_q_r, res_b_r, eC_r, cluster_idx=cidx_r, k_override=q+1),invG_p_r))
        V_T_cl = math.factorial(deriv)**2 * (V_T_cl_l + V_T_cl_r)[deriv, deriv]
        V_T_rb = math.factorial(deriv)**2 * (V_T_rb_l + V_T_rb_r)[deriv, deriv]
        se_tau_T_cl = np.sqrt(V_T_cl)
        se_tau_T_rb = np.sqrt(V_T_rb)

    V_tau_cl = scalepar**2*math.factorial(deriv)**2*(V_Y_cl_l+V_Y_cl_r)[deriv,deriv]
    V_tau_rb = scalepar**2*math.factorial(deriv)**2*(V_Y_rb_l+V_Y_rb_r)[deriv,deriv]
    se_tau_cl = np.sqrt(V_tau_cl)
    se_tau_rb = np.sqrt(V_tau_rb)

    tau = np.array([tau_cl, tau_bc, tau_bc]).reshape(-1,1)
    se  = np.array([se_tau_cl,se_tau_cl,se_tau_rb]).reshape(-1,1)
    z   =  tau/se
    pv  = 2*sct.norm.cdf(-np.abs(z))
    ci = np.column_stack((tau - quant*se,tau + quant*se))

    # Model description
    covs_label = "Covariate-adjusted " if covs is not None else ""
    n_clust = None
    n_clust_l = None
    n_clust_r = None
    if cluster is not None:
        n_clust_l = len(np.unique(eC_l.reshape(-1)))
        n_clust_r = len(np.unique(eC_r.reshape(-1)))
        n_clust = len(np.unique(np.concatenate([eC_l.reshape(-1), eC_r.reshape(-1)])))
    if fuzzy is None:
        if deriv==0: rdmodel = covs_label + "Sharp RD estimates using local polynomial regression."
        elif deriv==1: rdmodel = covs_label + "Sharp Kink RD estimates using local polynomial regression."
    else:
        if deriv==0: rdmodel = covs_label + "Fuzzy RD estimates using local polynomial regression."
        elif deriv==1: rdmodel = covs_label + "Fuzzy Kink RD estimates using local polynomial regression."
    if cluster is not None:
        rdmodel += " Std. errors are clustered (" + str(n_clust) + " clusters)."

    gamma_p_out = gamma_p if covs is not None else None

    # Preparing Output
    Estimate = pd.DataFrame(np.array([tau_cl,tau_bc, se_tau_cl, se_tau_rb]).reshape(1,-1),
                            columns = ["tau.us","tau.bc","se.us","se.rb"],
                            index = pd.Index(["Estimate"]))
    bws = pd.DataFrame(np.array([[h_l,h_r],[b_l,b_r]]),
                            columns = ["left","right"],
                            index = pd.Index(["h","b"]))
    label = ["Conventional","Bias-Corrected","Robust"]
    coef = pd.DataFrame(tau,
                        columns = ["Coeff"],
                        index = pd.Index(label))
    se = pd.DataFrame(se,
                      columns = ["Std. Err."],
                      index = pd.Index(label))
    z = pd.DataFrame(z,
                     columns = ["z-stat."],
                     index = pd.Index(label))
    pv = pd.DataFrame(pv,
                     columns = ["P>|z|"],
                     index = pd.Index(label))
    ci = pd.DataFrame(ci,
                      columns = ["CI Lower","CI Upper"],
                      index = pd.Index(["Conventional","Bias-Corrected","Robust"]))

    # First-stage (T) results for fuzzy designs; None for sharp.
    if fuzzy is not None:
        tau_T_vec = np.array([tau_T_cl, tau_T_bc, tau_T_bc]).reshape(-1, 1)
        se_T_vec  = np.array([se_tau_T_cl, se_tau_T_cl, se_tau_T_rb]).reshape(-1, 1)
        z_T_vec   = tau_T_vec / se_T_vec
        pv_T_vec  = 2 * sct.norm.cdf(-np.abs(z_T_vec))
        ci_T_arr  = np.column_stack((tau_T_vec - quant*se_T_vec,
                                     tau_T_vec + quant*se_T_vec))
        tau_T = pd.DataFrame(tau_T_vec, columns=["Coeff"], index=pd.Index(label))
        se_T  = pd.DataFrame(se_T_vec,  columns=["Std. Err."], index=pd.Index(label))
        z_T   = pd.DataFrame(z_T_vec,   columns=["z-stat."], index=pd.Index(label))
        pv_T  = pd.DataFrame(pv_T_vec,  columns=["P>|z|"], index=pd.Index(label))
        ci_T  = pd.DataFrame(ci_T_arr,  columns=["CI Lower","CI Upper"],
                             index=pd.Index(label))
    else:
        tau_T = se_T = z_T = pv_T = ci_T = None
        beta_T_p_l = beta_T_p_r = None

    return rdrobust_output(Estimate, bws, coef, se, z, pv, ci, beta_p_l, beta_p_r,
                    V_Y_cl_l, V_Y_cl_r, V_Y_rb_l, V_Y_rb_r,
                    [N_l,N_r], [N_h_l,N_h_r], [N_b_l,N_b_r], [M_l,M_r],
                    [tau_Y_cl_l,tau_Y_cl_r], [tau_Y_bc_l,tau_Y_bc_r],
                    c, p, q, [bias_l,bias_r], kernel_type, all,
                    vce_type, bwselect, level, masspoints,
                    rdmodel=rdmodel, n_clust=n_clust,
                    n_clust_l=n_clust_l, n_clust_r=n_clust_r, coef_covs=gamma_p_out,
                    tau_T=tau_T, se_T=se_T, z_T=z_T, pv_T=pv_T, ci_T=ci_T,
                    beta_T_p_l=beta_T_p_l, beta_T_p_r=beta_T_p_r)
