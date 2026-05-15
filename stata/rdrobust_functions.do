********************************************************************************
* RDROBUST STATA PACKAGE -- rdrobust_functions
* Authors: Sebastian Calonico, Matias D. Cattaneo, Max H. Farrell, Rocio Titiunik
********************************************************************************
*!version 11.0.0  2026-05-13

version 16.0

capture mata mata drop rdrobust_res()
mata
real matrix rdrobust_res(real matrix X, real matrix y, real matrix T, real matrix Z, real matrix m, real matrix hii, string vce, real scalar matches, dups, dupsid, real scalar d)
{
n = length(y)
dT=dZ=0
if (rows(T)>1) dT = 1
if (rows(Z)>1) dZ = cols(Z)
res = J(n,1+dT+dZ,.)		
if (vce=="nn") {
	for (pos=1; pos<=n; pos++) {
		rpos = dups[pos] - dupsid[pos]
		lpos = dupsid[pos] - 1
		while (lpos+rpos < min((matches,n-1))) {
			if (pos-lpos-1 <= 0) rpos = rpos + dups[pos+rpos+1]
			else if (pos+rpos+1>n) lpos = lpos + dups[pos-lpos-1]
			else if ((X[pos]-X[pos-lpos-1]) > (X[pos+rpos+1]-X[pos])) rpos = rpos + dups[pos+rpos+1]
			else if ((X[pos]-X[pos-lpos-1]) < (X[pos+rpos+1]-X[pos])) lpos = lpos + dups[pos-lpos-1]
			else {
				rpos = rpos + dups[pos+rpos+1]
				lpos = lpos + dups[pos-lpos-1]
			}
		}
		ind_J = (pos-lpos)::(pos+rpos)
		y_J   = sum(y[ind_J])-y[pos]
		Ji = length(ind_J)-1
		res[pos,1] = sqrt(Ji/(Ji+1))*(y[pos] :- y_J/Ji)
		if (dT==1) {
				T_J = sum(T[ind_J])-T[pos]
				res[pos,2] = sqrt(Ji/(Ji+1))*(T[pos] :- T_J/Ji)
		}
		if (dZ>0) {
			for (i=1; i<=dZ; i++) {
				Z_J = sum(Z[ind_J,i])-Z[pos,i]
				res[pos,1+dT+i] = sqrt(Ji/(Ji+1))*(Z[pos,i] :- Z_J/Ji)
			}
		}
	}		
}
else if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
 	 if (vce=="hc0") w = 1
	 if (vce=="hc1") w = sqrt(n/(n-d))
	 if (vce=="hc2") w = sqrt(1:/(1:-hii))
	 if (vce=="hc3") w =      1:/(1:-hii)
	 res[,1] = w:*(y-m[,1])
	 if (dT==1) res[,2] = w:*(T-m[,2])
	 if (dZ>0) {
		for (i=1; i<=dZ; i++) {
			res[,1+dT+i] = w:*(Z[,i]-m[,1+dT+i])
		}
	}
}
return(res)
}
mata mosave rdrobust_res(), replace
end

*************************************************************************************************************************************************************
capture mata mata drop rdrobust_kweight()
mata
real matrix rdrobust_kweight(real matrix X, real scalar c, real scalar h, string kernel)
{
u = (X:-c)/h
	if (kernel=="epanechnikov" | kernel=="epa") {
		w = (0.75:*(1:-u:^2):*(abs(u):<=1))/h
	}
	else if (kernel=="uniform" | kernel=="uni") {
		w = (0.5:*(abs(u):<=1))/h
	}
	else {
		w = ((1:-abs(u)):*(abs(u):<=1))/h
	}	
return(w)	
}
mata mosave rdrobust_kweight(), replace
end

*************************************************************************************************************************************************************
capture mata mata drop rdrobust_bw()
mata
real matrix rdrobust_bw(real matrix Y, real matrix X, real matrix T, real matrix Z, real matrix C, real matrix W, real scalar c, real scalar o, real scalar nu, real scalar o_B, real scalar h_V, real scalar h_B, real scalar scale, string vce, real scalar nnmatch, string kernel, dups, dupsid, covs_drop_coll, | string scalar cr_method, transmorphic vcache)
{
	// cr_method: "" | "cr1" | "crv2" | "crv3". Cluster path only.
	// vcache: optional asarray("string") keyed by "o_nu". T1 (2026-05-12):
	// V-fit depends only on (o, nu) at fixed h_V. Callers pass a per-side
	// asarray to share results across the 6-18 pilot calls per rdrobust.
	real scalar has_vcache, used_cache
	string scalar vkey
	real colvector packed
	if (args() < 20) cr_method = ""
	has_vcache = (args() >= 21)
	used_cache = 0
	vkey = sprintf("%g_%g", o, nu)
	dT = dZ = dC = eC = indC = 0
	dW = length(W)
	if (has_vcache) {
		if (asarray_contains(vcache, vkey)) {
			packed = asarray(vcache, vkey)
			V_V    = packed[1]
			BConst = packed[2]
			n_s    = packed[3]
			s      = packed[4..(3 + n_s)]
			if (rows(T)>1) dT = 1
			if (rows(Z)>1) dZ = cols(Z)
			if (rows(C)>1) dC = 1
			used_cache = 1
		}
	}
	if (used_cache == 0) {
	w = rdrobust_kweight(X, c, h_V, kernel)
	if (dW>1) {
		w = W:*w
	}
	ind_V = selectindex(w:> 0); eY = Y[ind_V];eX = X[ind_V];eW = w[ind_V]
	n_V = length(ind_V)
	D_V = eY
	// Q1: Vandermonde via successive multiplication.
	R_V = J(n_V,o+1,1)
	if (o >= 1) {
		u_V_tmp = eX :- c
		for (j=2; j<=(o+1); j++) R_V[.,j] = R_V[.,j-1] :* u_V_tmp
	}
	invG_V = cholinv(quadcross(R_V,eW,R_V))
	e_v = J((o+1),1,0); e_v[nu+1]=1
	s = 1
	if (rows(T)>1) {
		dT = 1
		eT = T[ind_V]
		D_V = D_V,eT
	}
	if (rows(Z)>1) {
		dZ = cols(Z)
		colsZ = (2+dT)::(2+dT+dZ-1)
		eZ   = Z[ind_V,]
		D_V  = D_V,eZ
		U    = quadcross(R_V:*eW,D_V)
		ZWD  = quadcross(eZ,eW,D_V)
		UiGU = quadcross(U[,colsZ],invG_V*U) 
		ZWZ  = ZWD[,colsZ] - UiGU[,colsZ] 
		ZWY  = ZWD[,1::1+dT] - UiGU[,1::1+dT] 		
		if (covs_drop_coll==0) gamma = cholinv(ZWZ)*ZWY
		if (covs_drop_coll==1) gamma =  invsym(ZWZ)*ZWY
		if (covs_drop_coll==2) gamma =    pinv(ZWZ)*ZWY
		s = 1 \ -gamma[,1]
	}
	if (rows(C)>1) {
		dC = 1
		eC =  C[ind_V] 
		indC = order(eC,1) 
	}
	beta_V = invG_V*quadcross(R_V:*eW,D_V)	
	if (dZ==0 & dT==1) {	
	    tau_Y = factorial(nu)*beta_V[nu+1,1]
		tau_T = factorial(nu)*beta_V[nu+1,2]
		s = (1/tau_T \ -(tau_Y/tau_T^2))
	}
	if (dZ>0 & dT==1) {	
				s_T = (1 \ -gamma[,2])
                tau_Y = factorial(nu)*s'*  vec((beta_V[nu+1,1],beta_V[nu+1,colsZ]))
				tau_T = factorial(nu)*s_T'*vec((beta_V[nu+1,2],beta_V[nu+1,colsZ]))
				s = (1/tau_T \ -(tau_Y/tau_T^2) \ -(1/tau_T)*gamma[,1] + (tau_Y/tau_T^2)*gamma[,2])
	}	
	dups_V=dupsid_V=predicts_V=0
	if (vce=="nn") {
		dups_V   = dups[ind_V]
		dupsid_V = dupsid[ind_V]
	}
	if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
		predicts_V=R_V*beta_V
		if (vce=="hc2" | vce=="hc3") {
			
			hii = rowsum((R_V*invG_V):*(R_V:*eW))
			
		}
	}	
	res_V = rdrobust_res(eX, eY, eT, eZ, predicts_V, hii, vce, nnmatch, dups_V, dupsid_V, o+1)
	// For CR2/CR3 the hat-matrix adjustment uses R_V * sqrt(W); pass invG_V
	// and sqrtRX_V. For other paths these are unused and can be empty.
	sqrtRX_V = (cr_method=="crv2" | cr_method=="crv3") ? R_V:*sqrt(eW) : J(0,0,.)
	invG_V_c = (cr_method=="crv2" | cr_method=="crv3") ? invG_V : J(0,0,.)
	V_V = (invG_V*rdrobust_vce(dT+dZ, s, R_V:*eW, res_V, eC, indC, invG_V_c, sqrtRX_V, cr_method, 0)*invG_V)[nu+1,nu+1]
	v = quadcross(R_V:*eW,((eX:-c):/h_V):^(o+1))
	Hp = J(o+1, 1, 1)
	for (j=1; j<=(o+1); j++) Hp[j] = h_V^((j-1))
	BConst = (diag(Hp)*(invG_V*v))[nu+1]
	if (has_vcache) {
		packed = (V_V \ BConst \ length(s) \ vec(s))
		asarray(vcache, vkey, packed)
	}
	}

	w = rdrobust_kweight(X, c, h_B, kernel)
	if (dW>1) {
		w = W:*w
	}
	ind = selectindex(w:> 0) 
	n_B = length(ind)
	eY = Y[ind];eX = X[ind];eW = w[ind]
	D_B = eY
	// Q1: Vandermonde via successive multiplication.
	R_B = J(n_B,o_B+1,1)
	if (o_B >= 1) {
		u_B_tmp = eX :- c
		for (j=2; j<=(o_B+1); j++) R_B[.,j] = R_B[.,j-1] :* u_B_tmp
	}
	invG_B = cholinv(quadcross(R_B,eW,R_B))
	if (dT==1) {
		eT = T[ind]
		D_B = D_B,eT
	}
	if (dZ>0) {
		eZ = Z[ind,]
		D_B = D_B,eZ
	}
	if (dC==1) {
		eC=C[ind]
		indC = order(eC,1) 
	}	
	beta_B = invG_B*quadcross(R_B:*eW,D_B)	
	BWreg=0
	if (scale>0) {
		e_B = J((o_B+1),1,0); e_B[o+2]=1
		dups_B=dupsid_B=hii=predicts_B=0
		if (vce=="nn") {
			dups_B   = dups[ind]
			dupsid_B = dupsid[ind]
		}
		if (vce=="hc0" | vce=="hc1" | vce=="hc2" | vce=="hc3") {
			predicts_B=R_B*beta_B
			if (vce=="hc2" | vce=="hc3") {
				
				hii = rowsum((R_B*invG_B):*(R_B:*eW))
				
			}
		}	
		res_B = rdrobust_res(eX, eY, eT, eZ, predicts_B, hii, vce, nnmatch, dups_B, dupsid_B,o_B+1)
		sqrtRX_B = (cr_method=="crv2" | cr_method=="crv3") ? R_B:*sqrt(eW) : J(0,0,.)
		invG_B_c = (cr_method=="crv2" | cr_method=="crv3") ? invG_B : J(0,0,.)
		V_B = (invG_B*rdrobust_vce(dT+dZ, s, R_B:*eW, res_B, eC, indC, invG_B_c, sqrtRX_B, cr_method, 0)*invG_B)[o+2,o+2]
		BWreg = 3*BConst^2*V_B
	}
	B =  sqrt(2*(o+1-nu))*BConst*(s'*beta_B[o+2,]')
	V = (2*nu+1)*h_V^(2*nu+1)*V_V
	R = scale*(2*(o+1-nu))*BWreg
	rate = 1/(2*o+3)
	return(V,B,R,rate)	
}
mata mosave rdrobust_bw(), replace
end

****************************************************
capture mata mata drop rdrobust_vce()
mata
real matrix rdrobust_vce(real scalar d, real matrix s, real matrix RX, real matrix res, real matrix C, real matrix ind, real matrix invG, real matrix sqrtRX, string scalar cr_method, real scalar k_override)
{
	// k_override: when > 0, overrides the (n-1)/(n-k) df correction's k for
	// CR1. Used when RX is a "score-like" matrix (e.g. Q_q) whose ncol is
	// smaller than the effective number of regressors estimated. Without it,
	// CR1 at h!=b under cluster used k=p+1 from cols(Q_q), inconsistent with
	// the q-regression direct path at h=b which uses k=q+1. Pass 0 to use
	// cols(RX) as before.
	// Compute the meat M of the sandwich estimator. Non-cluster path
	// (length(C)<=1) applies HC weighting via the pre-weighted residuals passed
	// in res; cluster paths implement CR1 (df correction), CR2 (Bell-McCaffrey
	// half-inverse hat-matrix adjustment), or CR3 (Pustejovsky-Tipton
	// full-inverse hat-matrix adjustment, Woodbury-optimized to O(k^3) per
	// cluster). invG (k x k) and sqrtRX (n x k = R*sqrt(W)) are only used by
	// the CR2/CR3 branches; callers may pass empty matrices (J(0,0,.)) when
	// not applicable. cr_method selects the cluster variant.

	real scalar k, n, g, i, j, l
	real matrix M
	real scalar has_sqrtRX, has_invG
	real matrix C_o, RX_o, res_o, sRX_o, info
	real matrix sRX_g, RX_g, ri, e_g, L_g, u_g, M_g, GmL, G
	real matrix F_sq, V, C_half, tC_half, Ru_g, adj, score_g, score_l
	real rowvector lambda
	real colvector sigma2, c_coef, sv
	real scalar u_g_l
	real matrix rng_idx
	real scalar w

	k = cols(RX)
	M = J(k,k,0)
	n = length(C)

	if (n <= 1) {
		// non-cluster: residuals already carry HC weights (hc0/hc1/hc2/hc3)
		if (d==0){
			SS = res:^2
			M  = quadcross(RX,SS,RX)
		}
		else {
			// T4: Sigma_{i,j} s_i s_j r_i r_j = (Sigma_l s_l r_l)^2
			// One weighted quadcross instead of an (d+1)^2 Mata loop.
			// s is a column vector (from rdrobust.ado: 1 \ -gamma_p[,1]...).
			real colvector r_comb
			r_comb = res * s
			M = quadcross(RX, r_comb :* r_comb, RX)
		}
		return(M)
	}

	// ---- cluster path ----
	has_sqrtRX = (rows(sqrtRX) > 1)
	has_invG   = (rows(invG)   > 1)
	C_o   = C[ind]
	RX_o  = RX[ind,]
	res_o = res[ind,]
	info  = panelsetup(C_o,1)
	g     = rows(info)
	if (has_sqrtRX) sRX_o = sqrtRX[ind,]
	else            sRX_o = J(0,0,.)

	// ---- CR3 (Pustejovsky-Tipton) via Sherman-Morrison-Woodbury ----
	if (cr_method == "crv3" & has_invG) {
		G = invsym(invG)    // k x k Gram matrix = R'WR
		w = 1                // CRV3 is approximately unbiased
		if (d==0) {
			for (i=1; i<=g; i++) {
				if (has_sqrtRX) {
					sRX_g = panelsubmatrix(sRX_o, i, info)
					e_g   = panelsubmatrix(res_o, i, info)[,1]
					L_g   = quadcross(sRX_g, sRX_g)
					u_g   = quadcross(sRX_g, sRX_g[,1]:*e_g)
				}
				else {
					RX_g  = panelsubmatrix(RX_o, i, info)
					e_g   = panelsubmatrix(res_o, i, info)[,1]
					L_g   = quadcross(RX_g, RX_g)
					u_g   = quadcross(RX_g, e_g)
				}
				GmL = G - L_g
				if (rank(GmL) == k) M_g = invsym(GmL)
				else                M_g = J(k,k,0)      // fallback to CR1-style
				score_g = u_g + L_g * (M_g * u_g)
				M = M + score_g * score_g'
			}
		}
		else {
			for (i=1; i<=g; i++) {
				if (has_sqrtRX) {
					sRX_g = panelsubmatrix(sRX_o, i, info)
					L_g   = quadcross(sRX_g, sRX_g)
				}
				else {
					RX_g = panelsubmatrix(RX_o, i, info)
					L_g  = quadcross(RX_g, RX_g)
				}
				ri  = panelsubmatrix(res_o, i, info)
				GmL = G - L_g
				if (rank(GmL) == k) M_g = invsym(GmL)
				else                M_g = J(k,k,0)
				sv  = J(k,1,0)
				for (l=1; l<=1+d; l++) {
					if (has_sqrtRX) u_g = quadcross(sRX_g, sRX_g[,1]:*ri[,l])
					else            u_g = quadcross(RX_g, ri[,l])
					score_l = u_g + L_g * (M_g * u_g)
					sv = sv + s[l] * score_l
				}
				M = M + sv * sv'
			}
		}
		return(w*M)
	}

	// ---- CR2 (Bell-McCaffrey) via eigendecomp of F_sq = C_half'*L_g*C_half ----
	if (cr_method == "crv2" & has_invG) {
		C_half  = cholesky(invG)    // lower-tri L: L*L' = invG
		tC_half = C_half'            // upper-tri L'
		w = 1                         // CRV2 is approximately unbiased
		if (d==0) {
			for (i=1; i<=g; i++) {
				if (has_sqrtRX) {
					sRX_g = panelsubmatrix(sRX_o, i, info)
					e_g   = panelsubmatrix(res_o, i, info)[,1]
					L_g   = quadcross(sRX_g, sRX_g)
					u_g   = quadcross(sRX_g, sRX_g[,1]:*e_g)
				}
				else {
					RX_g = panelsubmatrix(RX_o, i, info)
					e_g  = panelsubmatrix(res_o, i, info)[,1]
					L_g  = quadcross(RX_g, RX_g)
					u_g  = quadcross(RX_g, e_g)
				}
				F_sq = tC_half * L_g * C_half
				symeigensystem(F_sq, V=., lambda=.)
				sigma2 = lambda'
				for (j=1; j<=length(sigma2); j++) if (sigma2[j] < 0) sigma2[j] = 0
				c_coef = J(length(sigma2),1,0)
				for (j=1; j<=length(sigma2); j++) {
					if (sigma2[j] < 1e-14) c_coef[j] = 0
					else c_coef[j] = (1/sqrt(max((1-sigma2[j], 1e-8))) - 1) / sigma2[j]
				}
				Ru_g    = tC_half * u_g
				adj     = V * (c_coef :* (V' * Ru_g))
				score_g = u_g + L_g * (C_half * adj)
				M = M + score_g * score_g'
			}
		}
		else {
			for (i=1; i<=g; i++) {
				if (has_sqrtRX) {
					sRX_g = panelsubmatrix(sRX_o, i, info)
					L_g   = quadcross(sRX_g, sRX_g)
				}
				else {
					RX_g = panelsubmatrix(RX_o, i, info)
					L_g  = quadcross(RX_g, RX_g)
				}
				ri   = panelsubmatrix(res_o, i, info)
				F_sq = tC_half * L_g * C_half
				symeigensystem(F_sq, V=., lambda=.)
				sigma2 = lambda'
				for (j=1; j<=length(sigma2); j++) if (sigma2[j] < 0) sigma2[j] = 0
				c_coef = J(length(sigma2),1,0)
				for (j=1; j<=length(sigma2); j++) {
					if (sigma2[j] < 1e-14) c_coef[j] = 0
					else c_coef[j] = (1/sqrt(max((1-sigma2[j], 1e-8))) - 1) / sigma2[j]
				}
				sv = J(k,1,0)
				for (l=1; l<=1+d; l++) {
					if (has_sqrtRX) u_g = quadcross(sRX_g, sRX_g[,1]:*ri[,l])
					else            u_g = quadcross(RX_g, ri[,l])
					Ru_g    = tC_half * u_g
					adj     = V * (c_coef :* (V' * Ru_g))
					score_l = u_g + L_g * (C_half * adj)
					sv = sv + s[l] * score_l
				}
				M = M + sv * sv'
			}
		}
		return(w*M)
	}

	// ---- Default cluster path (CR1 with finite-sample df correction) ----
	// k_df: which k enters the (n-1)/(n-k) correction. Defaults to cols(RX)
	// unless caller passes a positive k_override.
	w = ((n-1) / (n - (k_override > 0 ? k_override : k))) * (g / (g-1))
	if (d==0){
		for (i=1; i<=g; i++) {
			RX_g = panelsubmatrix(RX_o, i, info)
			ri   = panelsubmatrix(res_o, i, info)
			u_g  = quadcross(RX_g, ri)'
			M = M + quadcross(u_g, u_g)
		}
	}
	else {
		for (i=1; i<=g; i++) {
			RX_g = panelsubmatrix(RX_o, i, info)
			ri   = panelsubmatrix(res_o, i, info)
			sv   = J(1+d, k, .)
			for (l=1; l<=1+d; l++) {
				sv[l,] = quadcross(RX_g, s[l]:*ri[,l])'
			}
			score_g = colsum(sv)
			M = M + quadcross(score_g, score_g)
		}
	}
	return(w*M)
}
mata mosave rdrobust_vce(), replace
end

capture mata mata drop rdrobust_vce_qq_cluster()
mata
real matrix rdrobust_vce_qq_cluster(
    real matrix Q, real matrix R_q, real colvector W_b, real matrix invG_q,
    real matrix res, real matrix C, real matrix ind,
    real scalar d, real matrix s, string scalar cr_method)
{
    // Decoupled cluster-robust variance for the bias-corrected estimator when h != b.
    // Sandwich uses Q (Q_q, mixed-bandwidth) but leverage comes from the q-regression
    // (bandwidth b). Matrix analog of using hii_q for HC2/HC3 in the same branch.
    // Returns the meat matrix M such that V_rb = invG_p * M * invG_p.
    //
    // Q:      n x k     Q_q_side matrix
    // R_q:    n x k_R   q-regression design (raw, no weights)
    // W_b:    n x 1     kernel weights at bandwidth b
    // invG_q: k_R x k_R q-regression inverse Gram
    // res:    n x (1+d) q-regression raw residuals
    // C:      n x 1     cluster id
    // ind:    n x 1     sort permutation of rows to match panelsetup
    // cr_method: "crv2" or "crv3"

    real scalar k, k_R, g, i, l, j, r
    real matrix M, C_o, R_o, res_o, Q_o, Wb_o, ri, info, G_q
    real matrix R_g, Q_g, sR_g, L_g, P_g
    real matrix V_sv, V_r, W_eig, T_mat, Lambda, GmL, M_g, A
    real colvector e_gl, Wb_g, u_Qq, u_q, sv, sigma2, sigma, tau, gamma, keep, sig, score_l
    real rowvector lambda, tau_row
    real colvector keep_idx

    k   = cols(Q)
    k_R = cols(R_q)
    M   = J(k, k, 0)

    C_o   = C[ind]
    R_o   = R_q[ind,]
    res_o = res[ind,]
    Q_o   = Q[ind,]
    Wb_o  = W_b[ind]
    info  = panelsetup(C_o, 1)
    g     = rows(info)
    G_q   = invsym(invG_q)

    for (i = 1; i <= g; i++) {
        Q_g  = panelsubmatrix(Q_o,   i, info)
        R_g  = panelsubmatrix(R_o,   i, info)
        Wb_g = panelsubmatrix(Wb_o,  i, info)
        sR_g = R_g :* sqrt(Wb_g)           // n_g x k_R
        L_g  = quadcross(sR_g, sR_g)       // k_R x k_R
        P_g  = quadcross(Q_g,  R_g)        // k x k_R
        ri   = panelsubmatrix(res_o, i, info)

        // T2: per-cluster Lambda (crv2) or M_g (crv3) depends only on
        // L_g / G_q / invG_q -- not on the residual column l. Compute
        // once per cluster, reuse across the l-loop (2026-05-12).
        real scalar use_simple
        use_simple = 0
        Lambda = J(0, 0, .)
        M_g    = J(0, 0, .)
        if (cr_method == "crv2") {
            symeigensystem(L_g, V_sv=., lambda=.)
            sigma2 = lambda'
            for (j = 1; j <= rows(sigma2); j++) if (sigma2[j] < 0) sigma2[j] = 0
            sigma = sqrt(sigma2)
            keep = (sigma :> max(sigma) * 1e-10)
            r = sum(keep)
            if (r == 0) {
                use_simple = 1
            }
            else {
                keep_idx = select((1..rows(sigma))', keep)
                V_r = V_sv[, keep_idx]
                sig = sigma[keep_idx]
                T_mat = diag(sig) * (V_r' * invG_q * V_r) * diag(sig)
                T_mat = 0.5 * (T_mat + T_mat')
                symeigensystem(T_mat, W_eig=., tau_row=.)
                tau = tau_row'
                for (j = 1; j <= rows(tau); j++) {
                    if (tau[j] < 0)       tau[j] = 0
                    if (tau[j] > 1-1e-10) tau[j] = 1 - 1e-10
                }
                gamma = 1 :/ sqrt(1 :- tau) :- 1
                A = (V_r :/ sig') * W_eig
                Lambda = A * (gamma :* A')
            }
        }
        else {
            GmL = G_q - L_g
            if (rank(GmL) == k_R) M_g = invsym(GmL)
            else                  M_g = J(k_R, k_R, 0)
        }

        sv = J(k, 1, 0)
        for (l = 1; l <= 1 + d; l++) {
            e_gl = ri[,l]
            u_Qq = quadcross(Q_g, e_gl)
            if (use_simple) {
                score_l = u_Qq
            }
            else if (cr_method == "crv2") {
                u_q     = quadcross(R_g, Wb_g :* e_gl)
                score_l = u_Qq + P_g * (Lambda * u_q)
            }
            else {
                u_q     = quadcross(R_g, Wb_g :* e_gl)
                score_l = u_Qq + P_g * (M_g * u_q)
            }
            sv = sv + s[l] * score_l
        }

        M = M + sv * sv'
    }

    return(M)
}
mata mosave rdrobust_vce_qq_cluster(), replace
end


capture mata mata drop rdrobust_groupid()
mata
real colvector rdrobust_groupid(real colvector x, real vector at)
{
        real scalar i, j
        real colvector result, p
        result = J(rows(x),1,.)
        j = length(at)
		
        for (i=rows(x); i>0; i--) {
                        if (x[i]>=.) continue
         
			for (; j>0; j--) {
				if (at[j]<=x[i]) break
			}
            if (j>0) result[i,1] = j
        }
		
        return(result)
}
mata mosave rdrobust_groupid(), replace 
end


capture mata mata drop rdrobust_median()
mata
real colvector rdrobust_median(real colvector x)
{	
  n = length(x)
  s = sort(x,1)
  if (mod(n,2)==1) {  	
	result = s[(n+1)/2]
  } 
  else{
    i = n/2
    med1 = s[i]
    med2 = s[i+1]
    result = (med1+med2)/2	
  }     
	return(result)  
}
mata mosave rdrobust_median(), replace 
end


capture mata mata drop rdrobust_collapse()
mata
real matrix rdrobust_collapse(real matrix x, real colvector id)
{	
	info  = panelsetup(id,1)
	g     = rows(info)
	mean = J(g,2,.)
	var  = J(g,1,.)
	
	for (i=1; i<=g; i++) {
			Xi       = panelsubmatrix(x,  i, info)
			mean[i,] = mean(Xi)
			var[i]   = variance(Xi[,2])
	}
		
		nobs =  (info[,2]-info[,1]) :+ 1
		result = nobs, mean, var	
		return(result)  		
}	
mata mosave rdrobust_collapse(), replace 
end
