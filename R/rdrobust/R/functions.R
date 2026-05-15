# Normalize the `covs` argument into a numeric vector/matrix.
#
# Accepted forms:
#   - one-sided formula (e.g. `~ z1 + z2 + I(z3^2)`):
#       processed with model.matrix; factors expand to contrasts; the
#       intercept column is dropped. Symbols not in `data` fall through
#       to `caller`.
#   - character vector of column names (requires `data = `):
#       selected as `data[, covs]` and coerced to a matrix.
#   - anything else (numeric vector, matrix, data.frame, NULL):
#       returned unchanged.
.rdrobust_resolve_covs <- function(covs, data, caller) {
  if (is.null(covs)) return(NULL)
  if (is.character(covs) && length(covs) == 0)
    stop("`covs` is empty. Supply at least one column name, or omit `covs`.",
         call. = FALSE)
  if ((is.matrix(covs) || is.data.frame(covs)) && ncol(covs) == 0)
    stop("`covs` has zero columns. Supply at least one covariate, or omit `covs`.",
         call. = FALSE)
  if (inherits(covs, "formula")) {
    environment(covs) <- caller
    mf <- stats::model.frame(covs,
                             data = if (!is.null(data)) data else caller,
                             na.action = stats::na.pass)
    out <- stats::model.matrix(covs, data = mf)
    if (ncol(out) >= 1 && colnames(out)[1] == "(Intercept)")
      out <- out[, -1, drop = FALSE]
    return(out)
  }
  if (is.character(covs)) {
    if (is.null(data))
      stop("Character `covs` requires `data = ` to resolve column names. ",
           "Pass a matrix or use `covs = df[, cols]` instead.", call. = FALSE)
    missing_cols <- setdiff(covs, colnames(data))
    if (length(missing_cols))
      stop("Column(s) not found in `data`: ",
           paste(sQuote(missing_cols), collapse = ", "), call. = FALSE)
    return(data.matrix(data[, covs, drop = FALSE]))
  }
  if (!(is.numeric(covs) || is.data.frame(covs)))
    stop("'covs' must be a one-sided formula, character vector of column names, ",
         "numeric matrix/vector, or data frame.", call. = FALSE)
  covs
}

# Shared input-validation helpers used by rdrobust / rdbwselect / rdplot.
# Each entry-point validates auxiliary-vector lengths against length(x)
# BEFORE subset filtering, so wrong-length inputs error explicitly instead
# of being silently recycled by `[subset]` / `[na.ok]` indexing.
.rdrobust_check_length <- function(arg, name, n) {
  if (is.null(arg)) return(invisible())
  m <- if (is.matrix(arg) || is.data.frame(arg)) nrow(arg) else length(arg)
  if (m != n)
    stop(sprintf("'%s' must have %s equal to length(x) (got %d, expected %d).",
                 name,
                 if (is.matrix(arg) || is.data.frame(arg)) "nrow" else "length",
                 m, n),
         call. = FALSE)
}

.rdrobust_check_subset <- function(subset, n) {
  if (is.null(subset)) return(invisible())
  if (is.logical(subset)) {
    if (length(subset) != n)
      stop(sprintf("Logical 'subset' must have length equal to length(x) (got %d, expected %d).",
                   length(subset), n),
           call. = FALSE)
  } else if (is.numeric(subset)) {
    if (any(!is.finite(subset)) || any(subset < 1) || any(subset > n) || any(subset != round(subset)))
      stop(sprintf("Numeric 'subset' must contain integer indices in 1..%d.", n),
           call. = FALSE)
  } else {
    stop("'subset' must be logical or integer.", call. = FALSE)
  }
}


# Fast equivalent of `split(seq_along(C), C)` for an atomic cluster vector.
# R's `split.default` calls `as.factor(C)` every invocation (~12ms for a
# 50k integer vector), which dominates the CR1/CR2/CR3 variance loop across
# the ~14 BW-pilot calls per rdrobust. This avoids that by sorting once
# and slicing on consecutive-equal runs; works for integer, numeric, and
# character clusters alike. Returns an UNNAMED list -- rdrobust_vce only
# iterates by index, so the level names are not needed.
.rdrobust_vander <- function(u, p) {
  # Q1: build the Vandermonde matrix [1, u, u^2, ..., u^p] by successive
  # multiplication. 3-5x faster than outer(u, 0:p, `^`) for small p; same
  # allocation cost. Called from rdrobust_bw + rdrobust.R hot paths.
  n <- length(u)
  if (p < 1) return(matrix(1, n, 1))
  out <- matrix(1, n, p + 1)
  for (j in 2:(p + 1)) out[, j] <- out[, j - 1] * u
  out
}

.rdrobust_cluster_idx <- function(C) {
  if (is.factor(C)) return(split(seq_along(C), C))
  n <- length(C)
  if (n == 0) return(list())
  ord <- order(C)
  C_s <- C[ord]
  chg <- c(TRUE, C_s[-1L] != C_s[-n])
  brk <- c(which(chg), n + 1L)
  g <- length(brk) - 1L
  out <- vector("list", g)
  for (i in seq_len(g)) out[[i]] <- ord[brk[i]:(brk[i + 1L] - 1L)]
  out
}


# Resolve a set of bare-name args against `data`, returning a named list.
# Used by rdrobust / rdbwselect / rdplot to enable `y = vote, x = margin,
# data = df` style calls. `mc` and `caller` must be supplied by the caller
# (they're frame-sensitive). Only the args actually present in `mc` are
# resolved; the rest are returned as NULL so callers can keep their
# default. The function deliberately does NOT touch the entry-point
# argument promises, so `y = vote` works even when `vote` is a column of
# `data` that doesn't exist in the caller environment.
.rdrobust_resolve_data <- function(mc, data, caller, args) {
  out <- setNames(vector("list", length(args)), args)
  for (a in args) {
    if (a %in% names(mc) && !is.null(mc[[a]]))
      out[[a]] <- eval(mc[[a]], envir = data, enclos = caller)
  }
  out
}

qrXXinv = function(x, ...) {
  G <- crossprod(x)
  R <- try(chol(G), silent = TRUE)
  if (inherits(R, "try-error")) ginv(G) else chol2inv(R)
}

rdrobust_kweight = function(X, c,  h,  kernel){
  u = (X-c)/h
  if (kernel %in% c("epanechnikov", "epa")) {
    w = (0.75*(1-u^2)*(abs(u)<=1))/h
  } else if (kernel %in% c("uniform", "uni")) {
    w = (0.5*(abs(u)<=1))/h
  } else {
    w = ((1-abs(u))*(abs(u)<=1))/h
  }
  return(w)	
}

rdrobust_res = function(X, y, T, Z, m, hii, vce, matches, dups, dupsid, d, crv3=FALSE, crv2=FALSE, has_cluster=FALSE) {
  n = length(y)
  dT=dZ=0
  if (!is.null(T)) dT = 1
  if (!is.null(Z)) dZ = ncol(Z)
  res = matrix(NA,n,1+dT+dZ)

  if (vce=="nn") {
    for (pos in 1:n) {
      rpos = dups[pos] - dupsid[pos]
      lpos = dupsid[pos] - 1
      while (lpos+rpos < min(c(matches,n-1))) {
        if (pos-lpos-1 <= 0) rpos = rpos + dups[pos+rpos+1]
        else if (pos+rpos+1>n) lpos = lpos + dups[pos-lpos-1]
        else if ((X[pos]-X[pos-lpos-1]) > (X[pos+rpos+1]-X[pos])) rpos = rpos + dups[pos+rpos+1]
        else if ((X[pos]-X[pos-lpos-1]) < (X[pos+rpos+1]-X[pos])) lpos = lpos + dups[pos-lpos-1]
        else {
          rpos = rpos + dups[pos+rpos+1]
          lpos = lpos + dups[pos-lpos-1]
        }
      }
      ind_J = max(c(0,(pos-lpos))):min(c(n,(pos+rpos)))
      y_J   = sum(y[ind_J])-y[pos]
      Ji = length(ind_J)-1
      res[pos,1] = sqrt(Ji/(Ji+1))*(y[pos] - y_J/Ji)
      if (!is.null(T)) {
        T_J = sum(T[ind_J])-T[pos]
        res[pos,2] = sqrt(Ji/(Ji+1))*(T[pos] - T_J/Ji)
      }
      if (!is.null(Z)) {
        for (i in 1:dZ) {
          Z_J = sum(Z[ind_J,i])-Z[pos,i]
          res[pos,1+dT+i] = sqrt(Ji/(Ji+1))*(Z[pos,i] - Z_J/Ji)
        }
      }
    }
  }
  else if (crv3 | crv2) {
    # CRV2/CRV3 mode: return raw (unweighted) residuals.
    # Cluster-level hat-matrix adjustment is handled in rdrobust_vce.
    res[,1] = y - m[,1]
    if (dT==1) res[,2] = T - m[,2]
    if (dZ>0) {
      for (i in 1:dZ) {
        res[,1+dT+i] = Z[,i] - m[,1+dT+i]
      }
    }
  }
  else {
    if (vce=="hc0") w = 1
    else if (vce=="hc1") w = if (has_cluster) 1 else sqrt(n/(n-d))
    else if (vce=="hc2") w = sqrt(1/pmax(1-hii, 1e-8))
    else                 w =      1/pmax(1-hii, 1e-8)
    res[,1] = w*(y-m[,1])
    if (dT==1) res[,2] = w*(T-m[,2])
    if (dZ>0) res[,(2+dT):(1+dT+dZ)] = w*(Z-m[,(2+dT):(1+dT+dZ)])
  }
  return(res)
}


rdrobust_bw = function(Y, X, T, Z, C, W, c, o, nu, o_B, h_V, h_B, scale, vce, nnmatch, kernel, dups, dupsid, covs_drop_coll, ginv.tol, vcache = NULL){
  # T1 (2026-05-12): the V-fit depends only on (o, nu) when h_V is held
  # fixed across calls on the same side. Callers (rdbwselect.R and
  # rdrobust.R) pass a per-side environment via `vcache` to share results
  # across the 6-18 pilot calls per rdrobust invocation. The B-fit
  # below still runs each call (its h_B differs per stage).
  dT = dZ = dC = 0
  crv3 = (vce=="crv3") & !is.null(C)
  crv2 = (vce=="crv2") & !is.null(C)
  key <- paste0(o, "_", nu)
  hit <- !is.null(vcache) && exists(key, envir = vcache, inherits = FALSE)
  if (hit) {
    cached <- get(key, envir = vcache, inherits = FALSE)
    V_V    <- cached$V_V
    BConst <- cached$BConst
    s      <- cached$s
    if (!is.null(T)) dT <- 1
    if (!is.null(Z)) dZ <- ncol(Z)
  } else {
  w = rdrobust_kweight(X, c, h_V, kernel)
  if (!is.null(W)) w = W*w

  ind_V = w> 0; eY = Y[ind_V];eX = X[ind_V];eW = w[ind_V]
  n_V = sum(ind_V)
  D_V = eY
  R_V = .rdrobust_vander(as.numeric(eX - c), o)
  invG_V = qrXXinv(R_V*sqrt(eW))
  s = 1
  eT=eC=eZ=NULL
  if (!is.null(T)) {
    dT = 1
    eT = T[ind_V]
    D_V = cbind(D_V,eT)
  }
  if (!is.null(Z)) {
    dZ = ncol(Z)
    eZ = Z[ind_V,,drop=FALSE]
    D_V = cbind(D_V,eZ)
    U = crossprod(R_V*eW,D_V)
    ZWD  = crossprod(eZ*eW,D_V)
    colsZ = (2+dT):max(c(2+dT+dZ-1,(2+dT)))
    UiGU =  crossprod(matrix(U[,colsZ],nrow=o+1),invG_V%*%U)
    ZWZ = ZWD[,colsZ] - UiGU[,colsZ]
    ZWY = ZWD[,1:(1+dT)] - UiGU[,1:(1+dT)]
    if (covs_drop_coll==1) {
      gamma = ginv(ZWZ, tol=ginv.tol)%*%ZWY
      }
    else {
      gamma = chol2inv(chol(ZWZ))%*%ZWY
    }
    s = c(1 , -gamma[,1])
  }
  if (!is.null(C)) {
    dC = 1
    eC =  C[ind_V]
  }
  beta_V = invG_V%*%crossprod(R_V*eW,D_V)
  if (is.null(Z) & !is.null(T)) {
    tau_Y = c(factorial(nu)*beta_V[nu+1,1])
    tau_T = c(factorial(nu)*beta_V[nu+1,2])
    s = c(1/tau_T , -(tau_Y/tau_T^2))
  }
  if (!is.null(Z) & !is.null(T)) {
    s_T = c(1 , -gamma[,2])
    tau_Y = c(factorial(nu)*t(s)%*%  c(beta_V[nu+1,1],beta_V[nu+1,colsZ]))
    tau_T = c(factorial(nu)*t(s_T)%*%c(beta_V[nu+1,2],beta_V[nu+1,colsZ]))
    s = c(1/tau_T , -(tau_Y/tau_T^2) , -(1/tau_T)*gamma[,1] + (tau_Y/tau_T^2)*gamma[,2])
  }
  dups_V=dupsid_V=predicts_V=hii=0

  if (vce=="nn") {
    dups_V   = dups[ind_V]
    dupsid_V = dupsid[ind_V]
  }

  if ((vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") & !crv3) {
    predicts_V = R_V%*%beta_V
    if (vce=="hc2" | vce=="hc3") {
      hii = rowSums((R_V%*%invG_V)*(R_V*eW))
    }
  } else if (crv3 | crv2) {
    predicts_V = R_V%*%beta_V
  }

        sqrtRX_V  = if (crv3 | crv2) R_V * sqrt(eW) else NULL
        invG_V_c  = if (crv3 | crv2) invG_V else NULL
        res_V = rdrobust_res(eX, eY, eT, eZ, predicts_V, hii, vce, nnmatch, dups_V, dupsid_V, o+1, crv3=crv3, crv2=crv2, has_cluster=!is.null(C))
        aux = rdrobust_vce(dT+dZ, s, R_V*eW, res_V, eC, invG=invG_V_c, sqrtRX=sqrtRX_V, crv2=crv2)
        V_V = (invG_V%*%aux%*%invG_V)[nu+1,nu+1]
        v = crossprod(R_V*eW,((eX-c)/h_V)^(o+1))
        Hp = 0
        for (j in 1:(o+1)) Hp[j] = h_V^((j-1))
        BConst = (Hp*(invG_V%*%v))[nu+1]
        if (!is.null(vcache)) {
          assign(key, list(V_V=V_V, BConst=BConst, s=s),
                 envir = vcache, inherits = FALSE)
        }
  }
        
        w = rdrobust_kweight(X, c, h_B, kernel)
        if (!is.null(W)) w = W*w
        ind = w> 0 
        n_B = sum(ind)
        eY = Y[ind];eX = X[ind];eW = w[ind]
        D_B = eY
        R_B = .rdrobust_vander(as.numeric(eX - c), o_B)
        invG_B = qrXXinv(R_B*sqrt(eW))
        eT=eC=eZ=NULL
        if (!is.null(T)) {
          eT = T[ind]
          D_B = cbind(D_B,eT)
        }
        if (!is.null(Z)) {
          eZ = Z[ind,,drop=FALSE]
          D_B = cbind(D_B,eZ)
        }
        if (!is.null(C)) {
          eC=C[ind]
        }	
        beta_B = invG_B%*%crossprod(R_B*eW,D_B)	
        BWreg=0
        if (scale>0) {
        e_B = matrix(0,(o_B+1),1); e_B[o+2]=1
        dups_B=dupsid_B=hii=predicts_B=0
        if (vce=="nn") {
          dups_B   = dups[ind]
          dupsid_B = dupsid[ind]
        }
        if ((vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") & !crv3) {
          predicts_B = R_B%*%beta_B
          if (vce=="hc2" | vce=="hc3") {
            hii = rowSums((R_B%*%invG_B)*(R_B*eW))
          }
        } else if (crv3 | crv2) {
          predicts_B = R_B%*%beta_B
        }
        sqrtRX_B = if (crv3 | crv2) R_B * sqrt(eW) else NULL
        invG_B_c = if (crv3 | crv2) invG_B else NULL
        res_B = rdrobust_res(eX, eY, eT, eZ, predicts_B, hii, vce, nnmatch, dups_B, dupsid_B, o_B+1, crv3=crv3, crv2=crv2, has_cluster=!is.null(C))
        V_B = (invG_B%*%rdrobust_vce(dT+dZ, s, R_B*eW, res_B, eC, invG=invG_B_c, sqrtRX=sqrtRX_B, crv2=crv2)%*%invG_B)[o+2,o+2]
		BWreg = 3*BConst^2*V_B
	}
	B =  sqrt(2*(o+1-nu))*BConst%*%(t(s)%*%(beta_B[o+2,]))
	V = (2*nu+1)*h_V^(2*nu+1)*V_V
	R = scale*(2*(o+1-nu))*BWreg
	rate = 1/(2*o+3)
  output = list(V=V,B=B,R=R,rate=rate)
  return(output)
}

rdrobust_vce = function(d, s, RX, res, C, invG=NULL, sqrtRX=NULL, crv2=FALSE,
                        k_override=NULL, cluster_idx=NULL) {
  # k_override: when supplied, overrides the (n-1)/(n-k) df correction's k
  # for CR1. Used when RX is a "score-like" matrix (e.g. Q_q) whose ncol is
  # smaller than the effective number of regressors estimated. Without it,
  # CR1 at h!=b under cluster used k=p+1 from ncol(Q_q), inconsistent with
  # the q-regression direct path at h=b which uses k=q+1.
  # cluster_idx: optional precomputed output of .rdrobust_cluster_idx(C);
  # callers that invoke rdrobust_vce multiple times with the same C should
  # compute it once and pass it through (T3 optimization, 2026-05-12).
  # Note: k is the matrix dimension; k_df is what enters the df correction.
  k    = ncol(as.matrix(RX))
  k_df = if (is.null(k_override)) k else k_override
  M = matrix(0,k,k)
  crv3_mode = !is.null(invG) && !is.null(C) && !crv2
  crv2_mode = !is.null(invG) && !is.null(C) &&  crv2

  if (is.null(C)) {
    w = 1
    if (d==0){
      M  = crossprod(c(res)*RX)
    }
    else {
      # T4: Σᵢⱼ sᵢsⱼ rᵢrⱼ factors as (Σₗ sₗ rₗ)². One weighted
      # crossprod instead of an (d+1)² R loop. With r_comb = res %*% s,
      # crossprod(RX * r_comb) returns Σₙ r_comb[n]² · RX[n,]ᵀ RX[n,].
      r_comb <- as.vector(res %*% s)
      M = crossprod(c(r_comb) * RX)
    }
  }
  else if (crv3_mode) {
    # CRV3: cluster-level hat-matrix adjustment (Pustejovsky & Tipton 2018).
    #
    # By the Sherman-Morrison-Woodbury identity:
    #   A_g = (I_{n_g} - sqrtRX_g %*% invG %*% t(sqrtRX_g))^{-1}
    #       = I + sqrtRX_g %*% (G - L_g)^{-1} %*% t(sqrtRX_g)
    # where G = solve(invG) and L_g = crossprod(sqrtRX_g).
    #
    # The cluster score collapses to a k-vector without forming A_g:
    #   u_g     = t(sqrtRX_g) %*% sqrte_g          (k x 1)
    #   score_g = u_g + L_g %*% solve(G - L_g) %*% u_g
    #           = (I_k + L_g %*% M_g) %*% u_g       (k x 1)
    #
    # This is O(k^3) per cluster rather than O(n_g^3).
    #
    # For V_rb (sqrtRX=NULL, RX=Q_q used as approximate hat design):
    # same formula but without sqrt-weight scaling; fallback to M_g=0 (CR1)
    # when G - L_g is near-singular.
    if (is.null(cluster_idx)) cluster_idx = .rdrobust_cluster_idx(C)
    g = length(cluster_idx)
    w = 1  # CRV3 is approximately unbiased; no df correction needed
    has_sqrtRX = !is.null(sqrtRX)
    G = tryCatch(solve(invG), error=function(e) ginv(invG))  # k x k Gram matrix

    if (d==0) {
      for (i in seq_along(cluster_idx)) {
        ind   = cluster_idx[[i]]
        if (has_sqrtRX) {
          sRX_g = sqrtRX[ind,,drop=FALSE]
          sw_g  = sRX_g[,1,drop=FALSE]
          e_g   = res[ind,1,drop=FALSE]
          L_g   = crossprod(sRX_g)
          u_g   = crossprod(sRX_g, sw_g * e_g)
          M_g   = tryCatch(solve(G - L_g), error=function(e) ginv(G - L_g))
          score_g = u_g + L_g %*% (M_g %*% u_g)
        } else {
          RX_g  = RX[ind,,drop=FALSE]
          e_g   = res[ind,1,drop=FALSE]
          L_g   = crossprod(RX_g)
          M_g   = tryCatch(solve(G - L_g), error=function(e) matrix(0, k, k))
          u_g   = crossprod(RX_g, e_g)
          score_g = u_g + L_g %*% (M_g %*% u_g)
        }
        M = M + tcrossprod(score_g)
      }
    }
    else {
      for (i in seq_along(cluster_idx)) {
        ind = cluster_idx[[i]]
        ri  = res[ind,,drop=FALSE]
        sv  = numeric(k)
        if (has_sqrtRX) {
          sRX_g = sqrtRX[ind,,drop=FALSE]
          sw_g  = sRX_g[,1,drop=FALSE]
          L_g   = crossprod(sRX_g)
          M_g   = tryCatch(solve(G - L_g), error=function(e) ginv(G - L_g))
          for (l in 1:(1+d)) {
            u_g_l   = crossprod(sRX_g, sw_g * ri[,l])
            score_l = u_g_l + L_g %*% (M_g %*% u_g_l)
            sv      = sv + as.vector(score_l) * s[l]
          }
        } else {
          RX_g = RX[ind,,drop=FALSE]
          L_g  = crossprod(RX_g)
          M_g  = tryCatch(solve(G - L_g), error=function(e) matrix(0, k, k))
          for (l in 1:(1+d)) {
            u_g_l   = crossprod(RX_g, ri[,l])
            score_l = u_g_l + L_g %*% (M_g %*% u_g_l)
            sv      = sv + as.vector(score_l) * s[l]
          }
        }
        M = M + tcrossprod(matrix(sv, k, 1))
      }
    }
  }
  else if (crv2_mode) {
    # CRV2: cluster-level half-inverse hat-matrix adjustment (Bell & McCaffrey 2002).
    #
    # Uses A_g = (I_{n_g} - H_g)^{-1/2} where H_g = sqrtRX_g invG sqrtRX_g'.
    # Since H_g has rank <= k, we work entirely in k-dimensional space:
    #
    #   F_g'F_g = chol(invG) %*% L_g %*% t(chol(invG))   (k x k, symmetric)
    #
    # whose eigenvalues sigma2_j are the squared singular values of F_g
    # (= eigenvalues of H_g).  The per-eigenvalue CRV2 coefficient is:
    #
    #   c_j = [(1 - sigma2_j)^{-1/2} - 1] / sigma2_j
    #
    # and the cluster score collapses to:
    #
    #   score_g = u_g + L_g t(R) V diag(c) V' R u_g      (k x 1)
    #
    # where R = chol(invG), V = eigenvectors of F_g'F_g, u_g = X_g'W_g e_g.
    # This is O(k^3) per cluster - no n_g x k matrix needed.
    if (is.null(cluster_idx)) cluster_idx = .rdrobust_cluster_idx(C)
    g         = length(cluster_idx)
    w         = 1   # CRV2 is approximately unbiased; no df correction needed
    has_sqrtRX = !is.null(sqrtRX)
    C_half    = tryCatch(chol(invG), error=function(e) NULL)
    tC_half   = if (!is.null(C_half)) t(C_half) else NULL

    crv2_ok   = !is.null(C_half)

    if (d==0) {
      for (i in seq_along(cluster_idx)) {
        ind = cluster_idx[[i]]
        if (has_sqrtRX) {
          sRX_g = sqrtRX[ind,,drop=FALSE]
          sw_g  = sRX_g[,1,drop=FALSE]
          e_g   = res[ind,1,drop=FALSE]
          L_g   = crossprod(sRX_g)
          u_g   = crossprod(sRX_g, sw_g * e_g)
        } else {
          RX_g  = RX[ind,,drop=FALSE]
          e_g   = res[ind,1,drop=FALSE]
          L_g   = crossprod(RX_g)
          u_g   = crossprod(RX_g, e_g)
        }
        if (crv2_ok) {
          F_sq   = C_half %*% L_g %*% tC_half
          eig_g  = eigen(F_sq, symmetric=TRUE)
          sigma2 = pmax(eig_g$values, 0)
          V_g    = eig_g$vectors
          c_coef = ifelse(sigma2 < 1e-14, 0,
                          (1/sqrt(pmax(1-sigma2, 1e-8)) - 1) / sigma2)
          Ru_g    = C_half %*% u_g
          adj     = V_g %*% (c_coef * crossprod(V_g, Ru_g))
          score_g = u_g + L_g %*% (tC_half %*% adj)
        } else {
          score_g = u_g
        }
        M = M + tcrossprod(score_g)
      }
    }
    else {
      for (i in seq_along(cluster_idx)) {
        ind = cluster_idx[[i]]
        ri  = res[ind,,drop=FALSE]
        sv  = numeric(k)
        if (has_sqrtRX) {
          sRX_g = sqrtRX[ind,,drop=FALSE]
          sw_g  = sRX_g[,1,drop=FALSE]
          L_g   = crossprod(sRX_g)
        } else {
          RX_g  = RX[ind,,drop=FALSE]
          L_g   = crossprod(RX_g)
        }
        if (crv2_ok) {
          F_sq     = C_half %*% L_g %*% tC_half
          eig_g    = eigen(F_sq, symmetric=TRUE)
          sigma2   = pmax(eig_g$values, 0)
          V_g      = eig_g$vectors
          c_coef   = ifelse(sigma2 < 1e-14, 0,
                            (1/sqrt(pmax(1-sigma2, 1e-8)) - 1) / sigma2)
          LtC_half = L_g %*% tC_half
          for (l in 1:(1+d)) {
            u_g_l   = if (has_sqrtRX) crossprod(sRX_g, sw_g * ri[,l]) else crossprod(RX_g, ri[,l])
            Ru_g_l  = C_half %*% u_g_l
            adj_l   = V_g %*% (c_coef * crossprod(V_g, Ru_g_l))
            score_l = u_g_l + LtC_half %*% adj_l
            sv      = sv + as.vector(score_l) * s[l]
          }
        } else {
          for (l in 1:(1+d)) {
            u_g_l = if (has_sqrtRX) crossprod(sRX_g, sw_g * ri[,l]) else crossprod(RX_g, ri[,l])
            sv    = sv + as.vector(u_g_l) * s[l]
          }
        }
        M = M + tcrossprod(matrix(sv, k, 1))
      }
    }
  }
  else {
    # Standard cluster-robust (CR1) with small-sample df correction
    if (is.null(cluster_idx)) cluster_idx = .rdrobust_cluster_idx(C)
    g = length(cluster_idx)
    n = length(C)
    w = ((n-1)/(n-k_df))*(g/(g-1))
    if (d==0){
      for (i in seq_along(cluster_idx)) {
        ind = cluster_idx[[i]]
        Xi  = RX[ind,,drop=FALSE]
        ri  = res[ind,,drop=FALSE]
        sv  = as.vector(crossprod(Xi, ri))
        M   = M + tcrossprod(matrix(sv, k, 1))
      }
    }
    else {
      for (i in seq_along(cluster_idx)) {
        ind = cluster_idx[[i]]
        Xi  = RX[ind,,drop=FALSE]
        ri  = res[ind,,drop=FALSE]
        sv  = numeric(k)
        for (l in 1:(1+d)) sv = sv + as.vector(crossprod(Xi, s[l]*ri[,l]))
        M = M + tcrossprod(matrix(sv, k, 1))
      }
    }
  }
  return(w*M)
}

J.fun = function(B,V,n) {ceiling((((2*B)/V)*n)^(1/3))}


# Decoupled cluster-robust variance for the bias-corrected estimator when h != b.
# The sandwich uses Q_q (mixed-bandwidth score matrix) but the leverage comes
# from the q-regression (bandwidth b), whose cluster hat matrix has eigvals in
# [0, 1) — matrix analog of using hii_q for HC2/HC3.
#
# Returns the meat matrix M such that V_rb = invG_p · M · invG_p.
# Back-conversion from transformed (sqrt(W_b)) to original space leaves:
#   CRV3: e_adj_g = e_g + R_q_g · (G_q - L_g)^{-1} · u_q_g
#   CRV2: e_adj_g = e_g + R_q_g · Lambda_g · u_q_g
# where u_q_g = R_q_g' W_b e_g and L_g = R_q_g' W_b R_q_g.
.rdrobust_vce_qq_cluster <- function(Q, R_q, W_b, invG_q, res, C,
                                     crv2 = FALSE, s = 1, d = 0,
                                     cluster_idx = NULL) {
  Q      <- as.matrix(Q)
  R_q    <- as.matrix(R_q)
  res    <- as.matrix(res)
  k      <- ncol(Q)
  k_R    <- ncol(R_q)
  M      <- matrix(0, k, k)
  if (is.null(cluster_idx)) cluster_idx <- .rdrobust_cluster_idx(C)
  G_q    <- tryCatch(solve(invG_q), error = function(e) ginv(invG_q))

  for (g in seq_along(cluster_idx)) {
    idx <- cluster_idx[[g]]
    Q_g    <- Q[idx, , drop = FALSE]
    R_g    <- R_q[idx, , drop = FALSE]
    Wb_g   <- W_b[idx]
    sR_g   <- R_g * sqrt(Wb_g)                      # n_g × k_R, transformed design
    L_g    <- crossprod(sR_g)                       # k_R × k_R  = R_q_g' W_b R_q_g
    P_g    <- crossprod(Q_g, R_g)                   # k × k_R     = Q_q_g' R_q_g

    # T2: per-cluster Lambda (CRV2) or M_g (CRV3) depends only on
    # L_g / G_q / invG_q — not on the residual column `l`. Compute
    # once per cluster, reuse across the l-loop.
    Lambda <- NULL
    M_g    <- NULL
    use_simple <- FALSE
    if (crv2) {
      eig_L <- tryCatch(eigen(L_g, symmetric = TRUE),
                        error = function(e) NULL)
      if (is.null(eig_L)) {
        use_simple <- TRUE
      } else {
        sigma2 <- pmax(eig_L$values, 0)
        sigma  <- sqrt(sigma2)
        V_sv   <- eig_L$vectors
        keep   <- sigma > max(sigma) * 1e-10
        r      <- sum(keep)
        if (r == 0) {
          use_simple <- TRUE
        } else {
          V_r   <- V_sv[, keep, drop = FALSE]
          sig   <- sigma[keep]
          Dsig  <- diag(sig, nrow = r, ncol = r)
          T_mat <- Dsig %*% crossprod(V_r, invG_q %*% V_r) %*% Dsig
          T_mat <- 0.5 * (T_mat + t(T_mat))
          eigT  <- eigen(T_mat, symmetric = TRUE)
          tau   <- pmin(pmax(eigT$values, 0), 1 - 1e-10)
          W_g   <- eigT$vectors
          gamma <- 1/sqrt(1 - tau) - 1
          Dinv  <- diag(1/sig, nrow = r, ncol = r)
          A     <- V_r %*% Dinv %*% W_g
          Lambda <- A %*% (gamma * t(A))
        }
      }
    } else {
      GmL <- G_q - L_g
      M_g <- tryCatch(solve(GmL), error = function(e) matrix(0, k_R, k_R))
    }

    sv <- numeric(k)
    for (l in seq_len(1 + d)) {
      e_gl <- res[idx, l]
      u_Qq <- crossprod(Q_g, e_gl)               # k × 1
      if (use_simple) {
        score_l <- as.vector(u_Qq)
      } else {
        u_q     <- crossprod(R_g, Wb_g * e_gl)   # k_R × 1
        adj     <- if (crv2) Lambda %*% u_q else M_g %*% u_q
        score_l <- as.vector(u_Qq + P_g %*% adj)
      }
      sv <- sv + score_l * s[l]
    }
    M <- M + tcrossprod(matrix(sv, k, 1))
  }
  M
}





covs_drop_fun <- function(z) {
  z    <- as.matrix(z)
  qr_z <- qr(z, tol = 1e-7)
  keep <- sort(qr_z$pivot[seq_len(qr_z$rank)])
  list(covs = z[, keep, drop = FALSE], ncovs = qr_z$rank)
}
