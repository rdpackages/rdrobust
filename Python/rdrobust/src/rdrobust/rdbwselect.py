#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 18:57:15 2021

@author: rmasini
"""

import numpy as np
import pandas as pd
from rdrobust.funs import *

def rdbwselect(y, x, c = None, fuzzy = None, deriv = None, p = None, q = None,
               covs = None, covs_drop = True, kernel = "tri", weights = None,
               bwselect = "mserd", vce = "nn", cluster = None, nnmatch = 3,
               scaleregul = 1, sharpbw = False, all = None, subset = None,
               masspoints = "adjust", bwcheck = None, bwrestrict = True,
               stdvars = False, prchk = True, data = None):
        
    
    '''
     Implements bandwidth selectors for local polynomial Regression Discontinuity (RD) point estimators and inference procedures developed in Calonico, Cattaneo and Titiunik (2014a), Calonico, Cattaneo and Farrell (2018), Calonico, Cattaneo, Farrell and Titiunik (2019) and Calonico, Cattaneo and Farrell (2020).
    
    Companion commands are: rdrobust for point estimation and inference procedures, and rdplot for data-driven RD plots (see Calonico, Cattaneo and Titiunik (2015a) for details).
    
    A detailed introduction to this command is given in Calonico, Cattaneo and Titiunik (2015b) and Calonico, Cattaneo, Farrell and Titiunik (2019). A companion Stata package is described in Calonico, Cattaneo and Titiunik (2014b).
    
    For more details, and related Stata and R packages useful for analysis of RD designs, visit https://rdpackages.github.io/
    
   
    Parameters
    ----------
    y	
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
    specifies the bandwidth selection procedure to be used. Options are:
    
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
    
    Note: MSE = Mean Square Error; CER = Coverage Error Rate. Default is bwselect=mserd. For details on implementation see Calonico, Cattaneo and Titiunik (2014a), Calonico, Cattaneo and Farrell (2018), and Calonico, Cattaneo, Farrell and Titiunik (2017), and the companion software articles.
    
    vce
    specifies the procedure used to compute the variance-covariance matrix estimator. Options are:

    nn for heteroskedasticity-robust nearest neighbor variance estimator.

    hc0 for heteroskedasticity-robust plug-in residuals variance estimator without weights.

    hc1 for heteroskedasticity-robust plug-in residuals variance estimator with hc1 weights.

    hc2 for heteroskedasticity-robust plug-in residuals variance estimator with hc2 weights.

    hc3 for heteroskedasticity-robust plug-in residuals variance estimator with hc3 weights.

    cr1 for cluster-robust plug-in residuals variance estimator with degrees-of-freedom weights. Requires cluster.

    cr2 for the Bell-McCaffrey (2002) bias-reduced cluster-robust variance estimator (CRV2). Requires cluster.

    cr3 for the Pustejovsky-Tipton (2018) cluster-robust variance estimator (CRV3), approximately unbiased with few clusters. Requires cluster.

    The CR2/CR3 leverage correction applies to both the conventional and the robust bias-corrected standard errors, including when the point-estimation bandwidth h differs from the bias-correction bandwidth b; in that case the cluster leverage is computed from the bias (b) regression.

    Default is vce=nn without cluster, and vce=cr1 with cluster. Auto-switching rules match rdrobust: hc2+cluster→cr2, hc3+cluster→cr3, and cr*+no-cluster falls back to the matching hc* variant.

    cluster
    indicates the cluster ID variable used for cluster-robust variance estimation. See vce option for auto-switching rules when cluster is combined with a heteroskedasticity-robust estimator.
    
    nnmatch	
    to be combined with for vce=nn for heteroskedasticity-robust nearest neighbor variance estimator with nnmatch indicating the minimum number of neighbors to be used. Default is nnmatch=3
    
    scaleregul	
    specifies scaling factor for the regularization term added to the denominator of the bandwidth selectors. Setting scaleregul = 0 removes the regularization term from the bandwidth selectors; default is scaleregul = 1.
    
    sharpbw	
    option to perform fuzzy RD estimation using a bandwidth selection procedure for the sharp RD model. This option is automatically selected if there is perfect compliance at either side of the threshold.
    
    all	
    if specified, rdbwselect reports all available bandwidth selection procedures.
    
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
    
    prchk
    internal check function.

    data
    optional pandas DataFrame. When supplied, `y`, `x`, `covs`, `cluster`,
    `fuzzy`, `weights`, and `subset` may be passed as column-name strings
    (or a list of names for `covs`), and a one-sided formula is accepted
    for `covs`.

    Returns
    -------
    N
    vector with sample sizes to the left and to the right of the cutoff.

    N_h
    vector with effective sample sizes (under the selected bandwidth) to the left and to the right of the cutoff.

    M
    vector with the number of unique observations to the left and to the right of the cutoff (when `masspoints` is not `off`).

    c
    cutoff value.

    p
    order of the local-polynomial used to construct the point-estimator.

    q
    order of the local-polynomial used to construct the bias-correction estimator.

    bws
    matrix containing the estimated bandwidths for each selected procedure. Method names are available as `bws.index`.

    bwselect
    bandwidth selection procedure employed.

    kernel
    kernel function used to construct the local-polynomial estimator(s).

    vce, masspoints
    variance estimator and mass-points option used.
    
    
    References
    ----------
    Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2018. On the Effect of Bias Estimation on Coverage Accuracy in Nonparametric Inference. Journal of the American Statistical Association, 113(522): 767-779.
    
    Calonico, S., M. D. Cattaneo, and M. H. Farrell. 2020. Optimal Bandwidth Choice for Robust Bias Corrected Inference in Regression Discontinuity Designs. Econometrics Journal, 23(2): 192-210.
    
    Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2017. rdrobust: Software for Regression Discontinuity Designs. Stata Journal 17(2): 372-404.
    
    Calonico, S., M. D. Cattaneo, M. H. Farrell, and R. Titiunik. 2019. Regression Discontinuity Designs using Covariates. Review of Economics and Statistics, 101(3): 442-451.
    
    Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014a. Robust Nonparametric Confidence Intervals for Regression-Discontinuity Designs. Econometrica 82(6): 2295-2326.
    
    Calonico, S., M. D. Cattaneo, and R. Titiunik. 2014b. Robust Data-Driven Inference in the Regression-Discontinuity Design. Stata Journal 14(4): 909-946.
    
    Calonico, S., M. D. Cattaneo, and R. Titiunik. 2015a. Optimal Data-Driven Regression Discontinuity Plots. Journal of the American Statistical Association 110(512): 1753-1769.
    
    Calonico, S., M. D. Cattaneo, and R. Titiunik. 2015b. rdrobust: An R Package for Robust Nonparametric Inference in Regression-Discontinuity Designs. R Journal 7(1): 38-51.
    
    Cattaneo, M. D., B. Frandsen, and R. Titiunik. 2015. Randomization Inference in the Regression Discontinuity Design: An Application to the Study of Party Advantages in the U.S. Senate. Journal of Causal Inference 3(1): 1-24.
    
    See Also
    rdrobust, rdplot
    
    Example
    -------
    >>> x = numpy.random.uniform(low=-1, high=1, size=1000)
    >>> y = 5+3*x+2*(x>=0) + numpy.random.uniform(size = 1000)
    >>> rdbwselect(y,x)
    '''

    loc = resolve_data_args(data, y=y, x=x, cluster=cluster,
                            fuzzy=fuzzy, weights=weights, subset=subset)
    y, x      = loc['y'], loc['x']
    cluster   = loc['cluster']
    fuzzy     = loc['fuzzy']
    weights   = loc['weights']
    subset    = loc['subset']
    covs = resolve_covs(covs, data)

    if prchk:
        x = np.array(x).reshape(-1,1)
        y = np.array(y).reshape(-1,1)

        # Validate auxiliary-vector lengths against length(x) BEFORE subset filter.
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
            subset = np.array(subset)
            x = x[subset]
            y = y[subset]
            if len(x) == 0:
                raise Exception("'subset' removed all observations.")

        if all is None: all = False

        if c is None: c = 0
        if not (np.isscalar(c) and np.isreal(c) and np.isfinite(c)):
            raise Exception(f"Cutoff 'c' must be a single finite numeric value (received: {c}).")

        if bwcheck is not None:
            if not (np.isscalar(bwcheck) and np.isfinite(bwcheck) and bwcheck >= 1
                    and bwcheck == round(bwcheck)):
                raise Exception("bwcheck must be a single positive integer")

        if (masspoints is not None and masspoints is not False
                and masspoints not in ("check", "adjust", "off", "")):
            raise Exception("masspoints must be one of 'check', 'adjust', 'off', or False")
        
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
        
        na_ok = complete_cases(x) & complete_cases(y)
        
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

        # Internal mapping
        if vce=='cr1': vce = 'hc1'
        elif vce=='cr2': vce = 'crv2'
        elif vce=='cr3': vce = 'crv3'

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
            if subset is not None: fuzzy = fuzzy[subset]
            na_ok = na_ok & complete_cases(fuzzy)
      
        if weights is not None:
            weights = np.array(weights).reshape(-1,1)
            if subset is not None: weights = weights[subset]
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
  
    ### reescaling
      
    x_q25, x_q75 = quantile_type2(x, [0.25, 0.75])
    x_iq = x_q75 - x_q25
    BWp = min(np.std(x, ddof =1),x_iq/1.349)
    x_sd = y_sd = 1
    if stdvars:
        y_sd = np.std(y,ddof = 1)
        x_sd = np.std(x,ddof = 1)
        y = y/y_sd
        x = x/x_sd
        c = c/x_sd
        BWp = min(1,(x_iq/x_sd)/1.349)
	  
    ###############################################
    X_l = x[x<c]
    X_r = x[x>=c]
    x_l_min = np.min(X_l)
    x_r_max = np.max(X_r)
    range_l = np.abs(c-x_l_min)
    range_r = np.abs(c-x_r_max)
      
    Y_l = y[x<c]
    Y_r = y[x>=c]
    N_l = len(X_l)   
    N_r = len(X_r)
    x_min = np.min(x)
    x_max = np.max(x)
    N = N_r + N_l
    
    M_l = N_l
    M_r = N_r
    
    # X_uniq_l/r and M_l/M_r are needed by both the masspoints inspection
    # branch AND the bwcheck branch. Hoist so bwcheck works even when
    # masspoints='off'.
    if masspoints=="check" or masspoints=="adjust" or bwcheck is not None:
      X_uniq_l = np.sort(np.unique(X_l))[::-1]
      X_uniq_r = np.unique(X_r)
      M_l = len(X_uniq_l)
      M_r = len(X_uniq_r)
    if masspoints=="check" or masspoints=="adjust":
      M = M_l + M_r
      mass_l = 1-M_l/N_l
      mass_r = 1-M_r/N_r				
      if mass_l>=0.2 or mass_r>=0.2:
        print("Mass points detected in the running variable.")
        if masspoints=="check": print("Try using option masspoints=adjust.")
        if bwcheck is None and masspoints=="adjust": bwcheck = 10
        
    covs_drop_coll = dZ = 0
    if covs_drop: covs_drop_coll = 1 
    
    if prchk:
        # Check for COLLINEARITY
        if covs is not None:
            dZ = ncol(covs)
            if covs_drop: 
                covs_check = covs_drop_fun(covs)
                if ncol(covs_check) < dZ:
                    covs = covs_check
                    dZ = ncol(covs)
    
    if kernel=="epanechnikov" or kernel=="epa":
        kernel_type = "Epanechnikov"
        C_c = 2.34
    elif kernel=="uniform" or kernel=="uni":
        kernel_type = "Uniform"
        C_c = 1.843
    else:
        kernel_type = "Triangular"
        C_c = 2.576
      
    vce_type = "NN"
    if vce=="hc0": vce_type = "HC0"
    if vce=="hc1": vce_type = "HC1"
    if vce=="hc2": vce_type = "HC2"
    if vce=="hc3": vce_type = "HC3"
    if vce=="crv2": vce_type = "CR2"
    if vce=="crv3": vce_type = "CR3"
    if cluster is not None and vce_type in ("NN","HC1"):
        vce_type = "CR1"

    #***********************************************************************
      
    Z_l = Z_r = T_l = T_r = C_l = C_r = g_l = g_r = None
      
    dups_l = np.zeros(N_l).astype(int)
    dupsid_l = np.zeros(N_l).astype(int)
    dups_r = np.zeros(N_r).astype(int)
    dupsid_r = np.zeros(N_r).astype(int)
    if vce == "nn":
        aux_l  = pd.DataFrame({'nn_l': np.ones(N_l), 'X_l': X_l })
        dups_l   = aux_l.groupby('X_l')['nn_l'].transform('sum').values.astype(int)
        dupsid_l = aux_l.groupby('X_l')['nn_l'].transform('cumsum').values.astype(int)

        aux_r  = pd.DataFrame({'nn_r': np.ones(N_r), 'X_r': X_r })
        dups_r   = aux_r.groupby('X_r')['nn_r'].transform('sum').values.astype(int)
        dupsid_r = aux_r.groupby('X_r')['nn_r'].transform('cumsum').values.astype(int)

    if covs is not None:
        Z_l  = covs[(x<c).reshape(-1),:]
        Z_r  = covs[(x>=c).reshape(-1),:]
    
    perf_comp = False
    if fuzzy is not None:
        T_l  = fuzzy[x<c]
        T_r  = fuzzy[x>=c]
        # Reject fully degenerate first stage (no variation, no jump). One-sided
        # non-compliance falls through to the perf_comp branch.
        if (np.var(T_l) == 0 and np.var(T_r) == 0
                and abs(np.mean(T_l) - np.mean(T_r)) < np.sqrt(np.finfo(float).eps)):
            raise ValueError(
                "Fuzzy RD: first-stage variable has no variation and no jump "
                "at the cutoff. The fuzzy estimator is not identified."
            )
        if np.var(T_l)==0 or np.var(T_r)==0: perf_comp = True
        if perf_comp or sharpbw==True:
            T_l = T_r = None
       
    if cluster is not None:
        C_l  = cluster[x<c]
        C_r = cluster[x>=c]
        g_l = len(np.unique(C_l))	
        g_r = len(np.unique(C_r))
    
    fw_l = fw_r = 0 
    if weights is not None:
        fw_l = weights[x<c]  
        fw_r = weights[x>=c]
                                                                                 
    #***********************************************************************
    bws_df, selected_bwselect = _rdbwselect_compute(
        Y_l=Y_l, Y_r=Y_r, X_l=X_l, X_r=X_r,
        T_l=T_l, T_r=T_r, Z_l=Z_l, Z_r=Z_r, C_l=C_l, C_r=C_r,
        fw_l=fw_l, fw_r=fw_r,
        dups_l=dups_l, dups_r=dups_r, dupsid_l=dupsid_l, dupsid_r=dupsid_r,
        N_l=N_l, N_r=N_r, N=N, M_l=M_l, M_r=M_r,
        M=(M_l + M_r) if masspoints in ("check", "adjust") else N,
        X_uniq_l=X_uniq_l if masspoints in ("check","adjust") else None,
        X_uniq_r=X_uniq_r if masspoints in ("check","adjust") else None,
        x_min=x_min, x_max=x_max, range_l=range_l, range_r=range_r, x_sd=x_sd,
        c=c, p=p, q=q, deriv=deriv, kernel=kernel, C_c=C_c, BWp=BWp,
        bwselect=bwselect, scaleregul=scaleregul, vce=vce, nnmatch=nnmatch,
        covs_drop_coll=covs_drop_coll,
        bwcheck=bwcheck, bwrestrict=bwrestrict, masspoints=masspoints,
        g_l=g_l, g_r=g_r, cluster_present=(cluster is not None), all=all,
    )
    bws = bws_df

    # Effective sample sizes for selected bandwidth
    h_sel_l = bws.iloc[0, 0]
    h_sel_r = bws.iloc[0, 1]
    w_h_l = rdrobust_kweight(X_l, c, h_sel_l, kernel)
    w_h_r = rdrobust_kweight(X_r, c, h_sel_r, kernel)
    N_h_l = int(np.sum(w_h_l > 0))
    N_h_r = int(np.sum(w_h_r > 0))

    return rdbwselect_output(bws, bwselect, kernel_type, p, q, c,
                            [N_l,N_r], [N_h_l,N_h_r], [M_l,M_r], vce_type, masspoints)


def _rdbwselect_compute(
    *, Y_l, Y_r, X_l, X_r, T_l, T_r, Z_l, Z_r, C_l, C_r, fw_l, fw_r,
    dups_l, dups_r, dupsid_l, dupsid_r,
    N_l, N_r, N, M_l, M_r, M, X_uniq_l, X_uniq_r,
    x_min, x_max, range_l, range_r, x_sd,
    c, p, q, deriv, kernel, C_c, BWp,
    bwselect, scaleregul, vce, nnmatch, covs_drop_coll,
    bwcheck, bwrestrict, masspoints, g_l, g_r, cluster_present, all,
):
    """Internal: compute bandwidths given pre-prepared side-specific data.

    Returns (bws_df, bwselect_name). Method names appear in bws_df.index.
    Called by rdbwselect() (public) and rdrobust() (to skip redundant prep).
    """

    c_bw = C_c*BWp*N**(-1/5)
    if masspoints=="adjust": c_bw = C_c*BWp*M**(-1/5)
          
    if bwrestrict:
        bw_max_l = abs(c-x_min)
        bw_max_r = abs(c-x_max)
        bw_max = max(bw_max_l, bw_max_r)
        c_bw = min(c_bw, bw_max)
          
    if bwcheck is not None:
        bwcheck_l = min(bwcheck, M_l)
        bwcheck_r = min(bwcheck, M_r)
        bw_min_l = np.abs(X_uniq_l-c)[bwcheck_l-1] + 1e-8
        bw_min_r = np.abs(X_uniq_r-c)[bwcheck_r-1] + 1e-8
        c_bw = max(c_bw, bw_min_l, bw_min_r)

    # Per-side V-fit caches: rdrobust_bw's V-fit depends only on (o, nu)
    # when h_V=c_bw is constant across all calls on a given side.
    vcache_l = {}
    vcache_r = {}

    #*** Step 1: d_bw
    C_d_l = (rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q+1, q+1,
                          q+2, c_bw, range_l, 0, vce, nnmatch,
                          kernel, dups_l, dupsid_l, covs_drop_coll,
                          _vcache=vcache_l))
    C_d_r = (rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q+1, q+1,
                          q+2, c_bw, range_r, 0, vce, nnmatch,
                          kernel, dups_r, dupsid_r, covs_drop_coll,
                          _vcache=vcache_r))
    #*** TWO bw 
    if (bwselect=="msetwo" or  bwselect=="certwo" or bwselect=="msecomb2" or
        bwselect=="cercomb2"  or all):		
        d_bw_l = (C_d_l[0]/(C_d_l[1]**2))**C_d_l[3]
        d_bw_r = (C_d_r[0]/(C_d_r[1]**2))**C_d_r[3]
        if bwrestrict:
            d_bw_l = min(d_bw_l, bw_max_l)
            d_bw_r = min(d_bw_r, bw_max_r)
        if bwcheck is not None:
            d_bw_l = max(d_bw_l, bw_min_l)
            d_bw_r = max(d_bw_r, bw_min_r)
        C_b_l  = (rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q, p+1, q+1,
                             c_bw, d_bw_l, scaleregul, vce, nnmatch, kernel,
                             dups_l, dupsid_l, covs_drop_coll,
                             _vcache=vcache_l))
        C_b_r  = (rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q, p+1, q+1,
                             c_bw, d_bw_r, scaleregul, vce, nnmatch, kernel,
                             dups_r, dupsid_r, covs_drop_coll,
                             _vcache=vcache_r))
        b_bw_l = (C_b_l[0]/(C_b_l[1]**2 + scaleregul*C_b_l[2]))**C_b_l[3]
        b_bw_r = (C_b_r[0]/(C_b_r[1]**2 + scaleregul*C_b_r[2]))**C_b_r[3]
        if bwrestrict:
            b_bw_l = min(b_bw_l, bw_max_l)
            b_bw_r = min(b_bw_r, bw_max_r)
        C_h_l  = (rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, p, deriv,
                              q, c_bw, b_bw_l, scaleregul, vce, nnmatch,
                              kernel, dups_l, dupsid_l, covs_drop_coll,
                              _vcache=vcache_l))
        C_h_r  = (rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, p, deriv,
                              q, c_bw, b_bw_r, scaleregul, vce, nnmatch,
                              kernel, dups_r, dupsid_r, covs_drop_coll,
                              _vcache=vcache_r))
        h_bw_l = (C_h_l[0]/(C_h_l[1]**2 + scaleregul*C_h_l[2]))**C_h_l[3]
        h_bw_r = (C_h_r[0]/(C_h_r[1]**2 + scaleregul*C_h_r[2]))**C_h_r[3]
        if bwrestrict:
            h_bw_l = min(h_bw_l, bw_max_l)
            h_bw_r = min(h_bw_r, bw_max_r)  
        h_msetwo_l = x_sd*h_bw_l
        h_msetwo_r = x_sd*h_bw_r
        b_msetwo_l = x_sd*b_bw_l
        b_msetwo_r = x_sd*b_bw_r  

    #*** SUM
    if (bwselect=="msesum" or bwselect=="cersum" or  bwselect=="msecomb1" or
         bwselect=="msecomb2" or  bwselect=="cercomb1" or bwselect=="cercomb2"
         or all):
        d_bw_s = ((C_d_l[0] + C_d_r[0])/(C_d_r[1] + C_d_l[1])**2)**C_d_l[3]
        if bwrestrict: d_bw_s = min(d_bw_s, bw_max)
        if bwcheck is not None: d_bw_s  =  max(d_bw_s, bw_min_l, bw_min_r)
        C_b_l  = (rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q, p+1, q+1,
                              c_bw, d_bw_s, scaleregul, vce, nnmatch, kernel,
                              dups_l, dupsid_l, covs_drop_coll,
                              _vcache=vcache_l))
        C_b_r  = (rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q, p+1, q+1,
                              c_bw, d_bw_s, scaleregul, vce, nnmatch, kernel,
                              dups_r, dupsid_r, covs_drop_coll,
                              _vcache=vcache_r))
        b_bw_s = ((C_b_l[0] + C_b_r[0])/((C_b_r[1] + C_b_l[1])**2 + scaleregul*(C_b_r[2]+C_b_l[2])))**C_b_l[3]

        if bwrestrict: b_bw_s = min(b_bw_s, bw_max)

        C_h_l  = (rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, p, deriv, q,
                              c_bw, b_bw_s, scaleregul, vce, nnmatch, kernel,
                              dups_l, dupsid_l, covs_drop_coll,
                              _vcache=vcache_l))
        C_h_r  = (rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, p, deriv, q,
                              c_bw, b_bw_s, scaleregul, vce, nnmatch, kernel,
                              dups_r, dupsid_r, covs_drop_coll,
                              _vcache=vcache_r))
        h_bw_s = ((C_h_l[0] + C_h_r[0])/((C_h_r[1] + C_h_l[1])**2 + scaleregul*(C_h_r[2]+C_h_l[2])))**C_h_l[3]
        if bwrestrict: h_bw_s = min(h_bw_s, bw_max)
        h_msesum = x_sd*h_bw_s
        b_msesum = x_sd*b_bw_s

    #*** RD
    if (bwselect=="mserd" or  bwselect=="cerrd" or bwselect=="msecomb1" or 
        bwselect=="msecomb2" or bwselect=="cercomb1" or bwselect=="cercomb2" 
        or bwselect=="" or all):
        d_bw_d = ((C_d_l[0] + C_d_r[0])/(C_d_r[1] - C_d_l[1])**2)**C_d_l[3]

        if bwrestrict: d_bw_d = min(d_bw_d, bw_max)
        if bwcheck is not None: d_bw_d = max(d_bw_d, bw_min_l, bw_min_r)
        C_b_l  = (rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, q, p+1, q+1,
                              c_bw, d_bw_d, scaleregul, vce, nnmatch, kernel,
                              dups_l, dupsid_l, covs_drop_coll,
                              _vcache=vcache_l))
        C_b_r  = (rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, q, p+1, q+1,
                              c_bw, d_bw_d, scaleregul, vce, nnmatch, kernel,
                              dups_r, dupsid_r, covs_drop_coll,
                              _vcache=vcache_r))
        b_bw_d = ((C_b_l[0] + C_b_r[0])/((C_b_r[1] - C_b_l[1])**2 + scaleregul*(C_b_r[2] + C_b_l[2])))**C_b_l[3]
        if bwrestrict: b_bw_d = min(b_bw_d, bw_max)
        C_h_l  = (rdrobust_bw(Y_l, X_l, T_l, Z_l, C_l, fw_l, c, p, deriv, q,
                              c_bw, b_bw_d, scaleregul, vce, nnmatch, kernel,
                              dups_l, dupsid_l, covs_drop_coll,
                              _vcache=vcache_l))
        C_h_r  = (rdrobust_bw(Y_r, X_r, T_r, Z_r, C_r, fw_r, c, p, deriv, q,
                              c_bw, b_bw_d, scaleregul, vce, nnmatch, kernel,
                              dups_r, dupsid_r, covs_drop_coll,
                              _vcache=vcache_r))
        h_bw_d = ((C_h_l[0] + C_h_r[0])/((C_h_r[1] - C_h_l[1])**2 + scaleregul*(C_h_r[2] + C_h_l[2])))**C_h_l[3]
        if bwrestrict: h_bw_d = min(h_bw_d, bw_max)
        h_mserd = x_sd*h_bw_d
        b_mserd = x_sd*b_bw_d
 
    if bwselect=="msecomb1" or bwselect=="cercomb1" or all: 
        h_msecomb1 = min(h_mserd,h_msesum)
        b_msecomb1 = min(b_mserd,b_msesum)
    if bwselect=="msecomb2" or bwselect=="cercomb2" or  all:
        h_msecomb2_l = np.median(np.array([h_mserd,h_msesum,h_msetwo_l]))
        h_msecomb2_r = np.median(np.array([h_mserd,h_msesum,h_msetwo_r]))
        b_msecomb2_l = np.median(np.array([b_mserd,b_msesum,b_msetwo_l]))
        b_msecomb2_r = np.median(np.array([b_mserd,b_msesum,b_msetwo_r]))
    cer_h = N**(-(p/((3+p)*(3+2*p))))
    if cluster_present: cer_h = (g_l+g_r)**(-(p/((3+p)*(3+2*p))))
    cer_b = 1
    if bwselect=="cerrd" or all:
        h_cerrd = h_mserd*cer_h
        b_cerrd = b_mserd*cer_b
    if bwselect=="cersum" or all:
        h_cersum = h_msesum*cer_h
        b_cersum=  b_msesum*cer_b
    if bwselect=="certwo" or all:
        h_certwo_l = h_msetwo_l*cer_h
        h_certwo_r = h_msetwo_r*cer_h
        b_certwo_l = b_msetwo_l*cer_b
        b_certwo_r = b_msetwo_r*cer_b		
    if bwselect=="cercomb1" or all:
        h_cercomb1 = h_msecomb1*cer_h
        b_cercomb1 = b_msecomb1*cer_b		
    if bwselect=="cercomb2" or all:
        h_cercomb2_l = h_msecomb2_l*cer_h
        h_cercomb2_r = h_msecomb2_r*cer_h
        b_cercomb2_l = b_msecomb2_l*cer_b
        b_cercomb2_r = b_msecomb2_r*cer_b
  
    if not all:
        if bwselect=="mserd" or bwselect=="":
            bws = np.array([h_mserd,h_mserd,b_mserd,b_mserd])
        if bwselect=="msetwo":
               bws = np.array([h_msetwo_l,h_msetwo_r,b_msetwo_l,b_msetwo_r])
        if bwselect=="msesum":
               bws = np.array([h_msesum,h_msesum,b_msesum,b_msesum])
        if bwselect=="msecomb1":
             bws = np.array([h_msecomb1,h_msecomb1,b_msecomb1,b_msecomb1])
        if bwselect=="msecomb2":
             bws = np.array([h_msecomb2_l,h_msecomb2_r,b_msecomb2_l,b_msecomb2_r]) 
        if bwselect=="cerrd":
                bws = np.array([h_cerrd,h_cerrd,b_cerrd,b_cerrd])
        if bwselect=="certwo":
               bws = np.array([h_certwo_l,h_certwo_r,b_certwo_l,b_certwo_r])
        if bwselect=="cersum":
               bws = np.array([h_cersum,h_cersum,b_cersum,b_cersum])
        if bwselect=="cercomb1":
             bws = np.array([h_cercomb1,h_cercomb1,b_cercomb1,b_cercomb1])
        if bwselect=="cercomb2":
             bws = np.array([h_cercomb2_l,h_cercomb2_r,b_cercomb2_l,b_cercomb2_r])        
        bws = pd.DataFrame(bws.reshape(1,-1), 
                           columns = ["h (left)","h (right)","b (left)","b (right)"],
                           index = pd.Index([bwselect]))
    if all:
        bwselect = "All"
        bws = nanmat(10,4)
        bws[0,:] = np.array([h_mserd,      h_mserd,      b_mserd,      b_mserd])
        bws[1,:] = np.array([h_msetwo_l,   h_msetwo_r,   b_msetwo_l,   b_msetwo_r])
        bws[2,:] = np.array([h_msesum,     h_msesum,     b_msesum,     b_msesum])
        bws[3,:] = np.array([h_msecomb1,   h_msecomb1,   b_msecomb1,   b_msecomb1])
        bws[4,:] = np.array([h_msecomb2_l, h_msecomb2_r, b_msecomb2_l, b_msecomb2_r]) 
        bws[5,:] = np.array([h_cerrd,      h_cerrd,      b_cerrd,      b_cerrd])
        bws[6,:] = np.array([h_certwo_l,   h_certwo_r,   b_certwo_l,   b_certwo_r])
        bws[7,:] = np.array([h_cersum,     h_cersum,     b_cersum,     b_cersum])
        bws[8,:] = np.array([h_cercomb1,   h_cercomb1,   b_cercomb1,   b_cercomb1])
        bws[9,:] = np.array([h_cercomb2_l, h_cercomb2_r, b_cercomb2_l, b_cercomb2_r])
        bws = pd.DataFrame(bws,
                           columns = ["h (left)","h (right)","b (left)","b (right)"],
                           index = pd.Index(["mserd","msetwo","msesum","msecomb1","msecomb2",
                                             "cerrd","certwo","cersum","cercomb1","cercomb2"]))
    return bws, bwselect
